# -*- coding: utf-8 -*-
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.cuda.amp import autocast, GradScaler

class Trainer:
    """
    Trainer
    - 支持 nn.DataParallel（多卡）
    - AMP 兼容老版本 autocast
    - microbatch 梯度累积
    """
    def __init__(self, model, optimizer, device=None,
                 use_amp=True, amp_dtype=torch.float16,
                 microbatch_size=None,
                 input_type='img', vis_dir=None, threshold=0.5, eval_sample=False):
        self.model = model
        self.optimizer = optimizer
        self.device = device
        self.use_amp = bool(use_amp and (device is not None and 'cuda' in str(device)))
        self.amp_dtype = amp_dtype
        self.microbatch_size = microbatch_size
        self.input_type = input_type
        self.vis_dir = vis_dir
        self.threshold = threshold
        self.eval_sample = eval_sample

        self.scaler = GradScaler(enabled=self.use_amp)

        if vis_dir is not None and not os.path.exists(vis_dir):
            os.makedirs(vis_dir)

    # ---------- helpers ----------
    def _use_dataparallel(self) -> bool:
        return isinstance(self.model, nn.DataParallel)

    def _to_tensor(self, x):
        """
        统一把 numpy -> Tensor；
        - 若 DataParallel：保持 CPU Tensor，交由 DP scatter 到各卡
        - 若 单卡 CUDA：迁移到 self.device
        """
        if isinstance(x, torch.Tensor):
            t = x.detach().to(dtype=torch.float32).contiguous()
        else:
            t = torch.as_tensor(x, dtype=torch.float32).contiguous()

        if self._use_dataparallel():
            # 让 DataParallel 自动 scatter CPU -> 多卡
            return t.cpu()
        else:
            # 单卡：直接搬到指定 device
            return t.to(self.device)

    def _iter_micro(self, data, label):
        """把 (B, ...) 按 microbatch_size 切分做累积；None 则不切分。"""
        if self.microbatch_size is None:
            yield data, label
            return
        B = data.shape[0]
        mb = int(self.microbatch_size)
        for s in range(0, B, mb):
            e = min(s + mb, B)
            yield data[s:e], label[s:e]

    def _autocast_ctx(self):
        """兼容旧版 AMP：老版没有 device_type / dtype 参数。"""
        try:
            return autocast(dtype=self.amp_dtype, enabled=self.use_amp)
        except TypeError:
            return autocast(enabled=self.use_amp)

    # ---------- core api ----------
    # train_utils/trainer.py

    def compute_loss(self, pred, target, *, save_debug=False, debug_prefix="debug"):
        """
        pred/target: (B, M, 3)
        使用单位化后的 MSE（等价于 1-cos 的稳定形式）
        """
        # —— 关键：把 target 搬到 pred 的 device/dtype —— #
        if not isinstance(pred, torch.Tensor):
            pred = torch.as_tensor(pred)
        if not isinstance(target, torch.Tensor):
            target = torch.as_tensor(target)

        # pred 来自 DataParallel 的主卡（cuda:0）；让 target 对齐它
        target = target.to(device=pred.device, dtype=pred.dtype, non_blocking=True)
        pred   = pred.to(device=pred.device, dtype=pred.dtype, non_blocking=True)

        if pred.shape != target.shape or pred.shape[-1] != 3:
            raise ValueError(f"Shape mismatch: pred {pred.shape}, target {target.shape}")

        if torch.isnan(pred).any() or torch.isinf(pred).any() or \
        torch.isnan(target).any() or torch.isinf(target).any():
            if save_debug:
                torch.save({"pred": pred.detach().cpu(),
                            "target": target.detach().cpu()},
                        f"{debug_prefix}_dump.pt")
            raise ValueError("[Loss Debug] NaN/Inf detected")

        pred_n   = F.normalize(pred,   dim=-1, eps=1e-8)
        target_n = F.normalize(target, dim=-1, eps=1e-8)
        loss = F.mse_loss(pred_n, target_n)
        return loss


    def train_step(self, data, label):
        """
        data:  (B, M, N, 3)  - numpy / tensor(任意设备)
        label: (B, M, 3)
        用 forward(model(...)) 触发 DataParallel 多卡拆分
        """
        self.model.train()
        self.optimizer.zero_grad(set_to_none=True)

        data  = self._to_tensor(data)
        label = self._to_tensor(label)

        total_loss = 0.0
        n_micro = 0

        for d_mb, y_mb in self._iter_micro(data, label):
            with self._autocast_ctx():
                pred = self.model(d_mb)         # (B_mb, M, 3)
                loss = self.compute_loss(pred, y_mb)

            if self.use_amp:
                self.scaler.scale(loss).backward()
            else:
                loss.backward()

            total_loss += float(loss.detach().cpu())
            n_micro += 1

        if self.use_amp:
            self.scaler.step(self.optimizer)
            self.scaler.update()
        else:
            self.optimizer.step()

        return total_loss / max(1, n_micro)

    def eval_step(self, data, label):
        """评估不反传。"""
        self.model.eval()
        data  = self._to_tensor(data)
        label = self._to_tensor(label)
        with torch.no_grad():
            with self._autocast_ctx():
                pred = self.model(data)
                loss = self.compute_loss(pred, label)
        return loss

    def evaluate(self, val_loader, sampling_size):
        val_sum = 0.0
        count = 0
        for i in range(val_loader.length()):
            for d, y in val_loader.iter_batches(i, sampling_size):
                loss = self.eval_step(d, y)
                B_mb = d.shape[0]
                val_sum += float(loss.detach().cpu()) * B_mb
                count += B_mb
        return val_sum / max(1, count)