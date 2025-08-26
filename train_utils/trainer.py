# 顶部 import
import math
import torch
import torch.nn.functional as F
import numpy as np

class Trainer:
    def __init__(self, model, optimizer, logger=None, cfg=None, device='cuda'):
        self.model = model
        self.opt = optimizer
        self.logger = logger
        self.device = device
        # 超参（可从 cfg 注入）
        self.aux_face_weight = getattr(cfg, "aux_face_weight", 0.1)     # 建议 0.05~0.2
        self.aux_center_only = getattr(cfg, "aux_center_only", True)    # 只监督中心面
        self.patch_detach = getattr(cfg, "patch_detach", False)         # True=完全切断（阶段2）
        self.patch_grad_scale = getattr(cfg, "patch_grad_scale", 0.0)   # 阶段3用，如 0.1

    # --- 核心：计算合成损失 ---
    def compute_loss_two_head(self, out, target):
        """
        out: (n_hat, n1) or n_hat
        target: (B,M,3)
        return: total_loss, dict(metrics)
        """
        tgt_n = F.normalize(target.to(self.device), dim=-1, eps=1e-8)

        if isinstance(out, tuple):
            n_hat, n1 = out
            n_hat = F.normalize(n_hat, dim=-1, eps=1e-8)
            n1    = F.normalize(n1,    dim=-1, eps=1e-8)

            # 主损失（最终输出）
            loss_main = F.mse_loss(n_hat, tgt_n)

            # 辅助损失（一层初稿）
            if self.aux_center_only:
                loss_face = F.mse_loss(n1[:, :1, :], tgt_n[:, :1, :])
            else:
                loss_face = F.mse_loss(n1, tgt_n)

            loss = loss_main + self.aux_face_weight * loss_face
            metrics = {
                "loss_main": float(loss_main.detach().cpu()),
                "loss_face": float(loss_face.detach().cpu()),
                "loss": float(loss.detach().cpu()),
            }
            return loss, metrics

        # 仅主输出（兼容老路径）
        n_hat = F.normalize(out, dim=-1, eps=1e-8)
        loss = F.mse_loss(n_hat, tgt_n)
        return loss, {"loss": float(loss.detach().cpu())}

    # --- 角度指标（度） ---
    def _angle_deg(self, pred, target):
        pn = F.normalize(pred,   dim=-1, eps=1e-8)
        tn = F.normalize(target, dim=-1, eps=1e-8)
        cos = torch.sum(pn * tn, dim=-1).clamp(-1 + 1e-7, 1 - 1e-7)
        return torch.acos(cos) * (180.0 / math.pi)  # (B,M)

    # --- 训练一步（支持微批循环的话，按你原来逻辑套进去即可） ---
    def train_step(self, data, target):
        self.model.train()
        if isinstance(data, np.ndarray):
            data = torch.from_numpy(data)
        if isinstance(target, np.ndarray):
            target = torch.from_numpy(target)
        data = data.to(self.device, non_blocking=True)
        target = target.to(self.device, non_blocking=True)

        self.opt.zero_grad(set_to_none=True)
        # 关键：在这里控制二层梯度回流
        out = self.model(
            data,
            patch_grad_s=self.patch_grad_scale,  # 例如 0.1 = 回流 10% 梯度
            detach_patch=self.patch_detach,      # True 时强制切断
        )
        loss, metrics = self.compute_loss_two_head(out, target)
        loss.backward()
        self.opt.step()

        # 可选日志
        if self.logger is not None and "loss_main" in metrics:
            self.logger.add_scalar("train/loss_main", metrics["loss_main"])
            self.logger.add_scalar("train/loss_face", metrics["loss_face"])
        return metrics["loss"]

    # --- 验证（只看主输出） ---
    def evaluate(self, dev_loader, sampling_size):
        self.model.eval()
        total_loss, total_cnt = 0.0, 0
        with torch.no_grad():
            for i in range(dev_loader.length()):
                data_batches, label_batches = dev_loader.generate_batch(i, sampling_size)
                for j in range(data_batches.shape[0]):
                    d = torch.from_numpy(data_batches[j]).to(self.device)
                    y = torch.from_numpy(label_batches[j]).to(self.device)
                    out = self.model(d, detach_patch=True)  # eval 时无梯度；detach=True 更保险
                    n_hat = out[0] if isinstance(out, tuple) else out
                    n_hat = F.normalize(n_hat, dim=-1, eps=1e-8)
                    y_n   = F.normalize(y,     dim=-1, eps=1e-8)
                    loss  = F.mse_loss(n_hat, y_n)
                    B = d.shape[0]
                    total_loss += float(loss.detach().cpu()) * B
                    total_cnt  += B
        return total_loss / max(1, total_cnt)

    # 可选：评估角度（主输出）
    def evaluate_angle(self, dev_loader, sampling_size, center_only=False):
        self.model.eval()
        total_ang, total_cnt = 0.0, 0
        with torch.no_grad():
            for i in range(dev_loader.length()):
                data_batches, label_batches = dev_loader.generate_batch(i, sampling_size)
                for j in range(data_batches.shape[0]):
                    d = torch.from_numpy(data_batches[j]).to(self.device)
                    y = torch.from_numpy(label_batches[j]).to(self.device)
                    out = self.model(d, detach_patch=True)
                    n_hat = out[0] if isinstance(out, tuple) else out
                    ang = self._angle_deg(n_hat, y)  # (B,M)
                    if center_only:
                        ang = ang[:, :1]
                    ang = ang.mean()
                    B = d.shape[0]
                    total_ang += float(ang.detach().cpu()) * B
                    total_cnt += B
        return total_ang / max(1, total_cnt)
