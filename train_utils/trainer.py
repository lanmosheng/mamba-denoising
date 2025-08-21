import torch
import torch.nn.functional as F
import os

class Trainer:
    def __init__(self, model, optimizer, device=None, input_type='img',
                 vis_dir=None, threshold=0.5, eval_sample=False):
        self.model = model
        self.optimizer = optimizer
        self.device = device
        self.input_type = input_type
        self.vis_dir = vis_dir
        self.threshold = threshold
        self.eval_sample = eval_sample
        if vis_dir is not None and not os.path.exists(vis_dir):
            os.makedirs(vis_dir)

    def train_step(self, data, label):
        self.model.train()
        self.optimizer.zero_grad()

        data  = data.to(self.device, dtype=torch.float32)
        label = label.to(self.device, dtype=torch.float32)

        # 两层 Mamba 前向（显式）
        face_features  = self.model.face_encoder(data)      # (B, M, d)
        patch_normals  = self.model.patch_encoder(face_features)  # (B, M, 3)

        loss = self.compute_loss(patch_normals, label)
        loss.backward()
        self.optimizer.step()
        return float(loss.detach().cpu())

    def eval_step(self, data, label):
        self.model.eval()
        data  = data.to(self.device, dtype=torch.float32)
        label = label.to(self.device, dtype=torch.float32)
        with torch.no_grad():
            face_features = self.model.face_encoder(data)
            patch_normals = self.model.patch_encoder(face_features)
            loss = self.compute_loss(patch_normals, label)
        return loss

    def compute_loss(self, pred, target, *, save_debug=True, debug_prefix="debug"):
        # 保证 tensor & 设备就绪（如果上游已是 tensor/同设备，这里基本是 no-op）
        pred   = pred   if isinstance(pred,   torch.Tensor) else torch.as_tensor(pred)
        target = target if isinstance(target, torch.Tensor) else torch.as_tensor(target)
        pred   = pred.to(self.device, dtype=torch.float32)
        target = target.to(self.device, dtype=torch.float32)

        # 形状/有效性检查
        if pred.shape != target.shape:
            raise ValueError(f"Shape mismatch: pred {pred.shape} vs target {target.shape}")
        if pred.shape[-1] != 3:
            raise ValueError(f"Expected last dim=3, got {pred.shape[-1]}")
        if torch.isnan(pred).any() or torch.isinf(pred).any() \
           or torch.isnan(target).any() or torch.isinf(target).any():
            if save_debug:
                torch.save({"pred": pred.detach().cpu(),
                            "target": target.detach().cpu()},
                           f"{debug_prefix}_dump.pt")
            raise ValueError("[Loss Debug] NaN/Inf detected in inputs")

        # 单位化后比较（注意 dim=-1）
        pred_n   = F.normalize(pred,   dim=-1, eps=1e-8)
        target_n = F.normalize(target, dim=-1, eps=1e-8)
        loss = F.mse_loss(pred_n, target_n)

        if torch.isnan(loss) or torch.isinf(loss):
            if save_debug:
                torch.save({"pred_n": pred_n.detach().cpu(),
                            "target_n": target_n.detach().cpu()},
                           f"{debug_prefix}_norm_dump.pt")
            raise ValueError("[Loss Debug] Loss is NaN/Inf")
        return loss
