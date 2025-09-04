# trainer_s1.py
import math
import torch
import torch.nn.functional as F

class TrainerS1:
    """
    仅用于第一层 Mamba（FaceEncoder）的训练与验证：
      - 训练：MSE（单位化后）
      - 验证：平均角度误差（度）
    依赖的数据加载器接口：
      loader.set_epoch(epoch)
      for batch in loader.iter_batches():
          batch['X']: (B, N, 3)  # 已做旋转规范化（不做 z-score）
          batch['Y']: (B, 3)     # 与 X 同一坐标系下的 GT（已旋转）
    """
    def __init__(self, model, optimizer, device="cuda", logger=None, grad_clip_norm=None, eps=1e-8):
        self.model = model
        self.opt = optimizer
        self.device = torch.device(device)
        self.logger = logger
        self.grad_clip_norm = grad_clip_norm
        self.eps = eps

    @staticmethod
    def _to_tensor(x, device):
        if isinstance(x, torch.Tensor):
            return x.to(device, non_blocking=True)
        import numpy as np
        if isinstance(x, np.ndarray):
            return torch.from_numpy(x).to(device, non_blocking=True)
        raise TypeError(f"Unsupported batch field type: {type(x)}")

    @staticmethod
    def _angle_deg(pred, target, eps=1e-8):
        # pred/target: (..., 3)
        pn = F.normalize(pred, dim=-1, eps=eps)
        tn = F.normalize(target, dim=-1, eps=eps)
        cos = torch.sum(pn * tn, dim=-1).clamp(-1 + 1e-7, 1 - 1e-7)
        return torch.acos(cos) * (180.0 / math.pi)  # (...,)

    def train_one_epoch(self, loader, epoch:int):
        """
        训练一轮；返回该轮的平均 MSE（单位化后）
        注意：会调用 loader.set_epoch(epoch) 以启用“跨 epoch 不重叠切片”的策略
        """
        self.model.train()
        loader.set_epoch(epoch)

        total_loss, total_cnt = 0.0, 0
        step = 0

        for batch in loader.iter_batches():
            X = self._to_tensor(batch['X'], self.device)  # (B,N,3)
            Y = self._to_tensor(batch['Y'], self.device)  # (B,3)

            self.opt.zero_grad(set_to_none=True)

            out = self.model(X)
            n_pred = out[0] if isinstance(out, tuple) else out  # (B,3)

            # 训练用 MSE（单位化后）
            n_pred = F.normalize(n_pred, dim=-1, eps=self.eps)
            Y_n    = F.normalize(Y,      dim=-1, eps=self.eps)
            loss = F.mse_loss(n_pred, Y_n)

            loss.backward()
            if self.grad_clip_norm is not None:
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.grad_clip_norm)
            self.opt.step()

            B = X.shape[0]
            total_loss += float(loss.detach().cpu()) * B
            total_cnt  += B
            step += 1

            if self.logger is not None:
                self.logger.add_scalar("train/mse", float(loss.detach().cpu()))

        return total_loss / max(1, total_cnt)

    def validate_angle(self, loader):
        """
        计算平均角度误差（度）。不修改 loader 的 epoch 切片策略：
        - 若想全量评估，请构造一个 slice_ratio=1.0 的验证 loader
        - 或者外部先 set_epoch(fixed) 来固定切片
        """
        self.model.eval()
        ang_sum, ang_cnt = 0.0, 0
        with torch.no_grad():
            for batch in loader.iter_batches():
                X = self._to_tensor(batch['X'], self.device)  # (B,N,3)
                Y = self._to_tensor(batch['Y'], self.device)  # (B,3)

                out = self.model(X, )
                n_pred = out[0] if isinstance(out, tuple) else out  # (B,3)

                # 平均角度误差（单位：度）
                ang = self._angle_deg(n_pred, Y, eps=self.eps)  # (B,)
                ang_sum += float(ang.sum().detach().cpu())
                ang_cnt += int(ang.numel())

        mean_deg = ang_sum / max(1, ang_cnt)
        if self.logger is not None:
            self.logger.add_scalar("val/angle_deg", mean_deg)
        return mean_deg
