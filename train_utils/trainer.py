# -*- coding: utf-8 -*-
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
        # 训练超参（与现有保持一致）
        self.aux_face_weight = getattr(cfg, "aux_face_weight", 0.1)
        self.patch_detach = getattr(cfg, "patch_detach", False)
        self.patch_grad_scale = getattr(cfg, "patch_grad_scale", 0.0)

        # 前处理配置（旋转 + L2 归一化）
        self.eps              = getattr(cfg, "eps", 1e-6) if cfg is not None else 1e-6

    # ---------------- 前处理：旋转(+X) + L2归一化 ----------------
    def _preprocess(self, X, Y):
        """
        输入:
          X: (B,M,N,3)  原始 LSD
          Y: (B,M,3)    原始 GT 法向（全局坐标）
        输出:
          X_rot: (B,M,N,3)  旋到 +X 的 LSD（送模型）
          Y_glb: (B,M,3)    保持“全局坐标”的 GT（不旋转，用于方案A聚合）
          R:     (B,3,3)    每个样本的旋转矩阵（a->+X）
          Y_rot: (B,M,3)    旋转后的 GT（评估/非聚合路径复用）
        """
        B, M, N, _ = X.shape
        centers = X[:, :, 0, :]  # (B,M,3)

        anchor = centers[:, 0, :]  # (B,3) patch中心面
        

        # 兜底：anchor 异常时改用有限均值；再不行用 [1,0,0]
        bad = (torch.isnan(anchor).any(dim=-1) | torch.isinf(anchor).any(dim=-1) |
                (torch.linalg.norm(anchor, dim=-1) < 1e-6))
        if bad.any():
            Xflat = X.view(B, -1, 3)
            mask  = torch.isfinite(Xflat).all(dim=-1, keepdim=True).float()
            mean  = (Xflat * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1.0)
            anchor = torch.where(bad.unsqueeze(-1), mean, anchor)

        ex = torch.tensor([1.0, 0.0, 0.0], device=X.device, dtype=X.dtype).expand(X.shape[0], -1)
        a = anchor / torch.linalg.norm(anchor, dim=-1, keepdim=True).clamp_min(self.eps)
        b = ex

        # Rodrigues 批量构造 R(a->b)
        v = torch.cross(a, b, dim=-1)           # (B,3)
        c = (a * b).sum(dim=-1)                 # (B,)
        s = torch.linalg.norm(v, dim=-1)        # (B,)
        I = torch.eye(3, device=X.device, dtype=X.dtype).unsqueeze(0).expand(B, -1, -1)

        # 平行/反向分支
        near0 = (s < 1e-12)
        # 反向：选与 a 正交的轴
        alt_x = torch.tensor([1.0, 0.0, 0.0], device=X.device, dtype=X.dtype).expand_as(a)
        alt_y = torch.tensor([0.0, 1.0, 0.0], device=X.device, dtype=X.dtype).expand_as(a)
        use_alt_y = (a[:, 0].abs() > 0.9).unsqueeze(-1)  # (B,1)
        alt = torch.where(use_alt_y, alt_y, alt_x)       # (B,3)
        axis_anti = torch.cross(a, alt, dim=-1)
        axis_anti = axis_anti / torch.linalg.norm(axis_anti, dim=-1, keepdim=True).clamp_min(self.eps)
        zeros = torch.zeros(B, device=X.device, dtype=X.dtype)
        K180 = torch.stack([
            torch.stack([zeros, -axis_anti[:, 2], axis_anti[:, 1]], dim=-1),
            torch.stack([axis_anti[:, 2], zeros, -axis_anti[:, 0]], dim=-1),
            torch.stack([-axis_anti[:, 1], axis_anti[:, 0], zeros], dim=-1)
        ], dim=1)
        R_anti = I + 2.0 * (K180 @ K180)

        # 一般情况
        axis = v / s.clamp_min(self.eps).unsqueeze(-1)
        K = torch.stack([
            torch.stack([zeros, -axis[:, 2], axis[:, 1]], dim=-1),
            torch.stack([axis[:, 2], zeros, -axis[:, 0]], dim=-1),
            torch.stack([-axis[:, 1], axis[:, 0], zeros], dim=-1)
        ], dim=1)
        sin_term = s.view(-1, 1, 1)
        one_minus_c = (1.0 - c).view(-1, 1, 1)
        R_general = I + K * sin_term + (K @ K) * one_minus_c

        c_pos = (c > 0)
        R = torch.where(near0.view(-1, 1, 1) & c_pos.view(-1, 1, 1), I, R_general)
        R = torch.where(near0.view(-1, 1, 1) & (~c_pos).view(-1, 1, 1), R_anti, R)

        # 应用旋转：X_rot = R^T @ X，Y_glb 不旋转（保留全局）
        XT = X.view(B, -1, 3)
        RT = R.transpose(1, 2)
        X_rot = torch.einsum('bij, bkj -> bki', RT, XT).view_as(X)
        Y_glb = Y
        # 提供一个 Y_rot（用于非聚合路径/评估）
        YT = Y.view(B, -1, 3)
        Y_rot = torch.einsum('bij, bkj -> bki', RT, YT).view_as(Y)

        # L2 归一化（按配置）
        Xn = torch.linalg.norm(X_rot, dim=-1, keepdim=True).clamp_min(self.eps)
        X_rot = X_rot / Xn
    
        Yn = torch.linalg.norm(Y_glb, dim=-1, keepdim=True).clamp_min(self.eps)
        Y_glb = Y_glb / Yn
        Yn_r = torch.linalg.norm(Y_rot, dim=-1, keepdim=True).clamp_min(self.eps)
        Y_rot = Y_rot / Yn_r

        return X_rot, Y_glb, R, Y_rot

    
    # ---------------- 方案A：反旋回 + 按面聚合（向量平均） ----------------
    @staticmethod
    def _unrotate_vec(pred_rot, R):
        """pred_rot: (B,M,3) in rotated coords; R: (B,3,3); return pred in GLOBAL coords."""
        return torch.einsum('bij, bmj -> bmi', R, pred_rot)

    @staticmethod
    def _agg_by_face(vec, face_idx):
        """
        vec: (B,M,3)  全局坐标下的向量
        face_idx: (B,M)  全局面号（Long）
        return:
          mean_vec: (U,3)  去重后的“按面平均”结果
        """
        B, M, _ = vec.shape
        flat_v = vec.reshape(B * M, 3)
        flat_i = face_idx.reshape(B * M)
        # 去重聚合
        uniq, inv = torch.unique(flat_i, return_inverse=True)
        U = uniq.shape[0]
        sums = torch.zeros(U, 3, device=vec.device, dtype=vec.dtype)
        cnts = torch.zeros(U, 1, device=vec.device, dtype=vec.dtype)
        sums.index_add_(0, inv, flat_v)
        cnts.index_add_(0, inv, torch.ones(B * M, 1, device=vec.device, dtype=vec.dtype))
        mean_vec = sums / cnts.clamp_min(1.0)
        return mean_vec

    def compute_loss_two_head_faceagg(self, out, target_glb, R, face_idx):
        """
        out: (n_hat_rot, n1_rot) or n_hat_rot   —— 模型在“旋转坐标系”的输出
        target_glb: (B,M,3) —— GT（全局坐标），来自 Loader 的原始 Y（已 L2）
        R: (B,3,3) —— 每个样本的旋转矩阵（a->+X）
        face_idx: (B,M) —— 全局面号（Long）
        """
        # 主输出：反旋回到全局系
        if isinstance(out, tuple):
            n_hat_rot, n1_rot = out
        else:
            n_hat_rot, n1_rot = out, None

        n_hat_glb = self._unrotate_vec(n_hat_rot, R)           # (B,M,3)
        # 按面聚合：对同一面在不同 patch 的预测向量取平均
        pred_face = self._agg_by_face(n_hat_glb, face_idx)     # (U,3)
        gt_face   = self._agg_by_face(target_glb, face_idx)    # (U,3)

        # 归一化后做 MSE
        pred_face = F.normalize(pred_face, dim=-1, eps=1e-8)
        gt_face   = F.normalize(gt_face,   dim=-1, eps=1e-8)
        loss_main = F.mse_loss(pred_face, gt_face)

        # 辅助头
        if n1_rot is not None and self.aux_face_weight > 0:
            # 辅助头默认保留“中心面监督”更稳（也可改成面级聚合）
            n1_glb = self._unrotate_vec(n1_rot, R)  # (B,M,3)
            
            # 面级聚合版本（可选）
            pred_aux = self._agg_by_face(n1_glb, face_idx)
            gt_aux   = self._agg_by_face(target_glb, face_idx)
            loss_face = F.mse_loss(
                F.normalize(pred_aux, dim=-1, eps=1e-8),
                F.normalize(gt_aux,   dim=-1, eps=1e-8),
            )
        else:
            loss_face = torch.tensor(0.0, device=loss_main.device, dtype=loss_main.dtype)

        loss = loss_main + self.aux_face_weight * loss_face
        metrics = {
            "loss_main": float(loss_main.detach().cpu()),
            "loss_face": float(loss_face.detach().cpu()),
            "loss": float(loss.detach().cpu()),
        }
        return loss, metrics

    # ---------------- 角度（度） ----------------
    def _angle_deg(self, pred, target):
        pn = F.normalize(pred,   dim=-1, eps=1e-8)
        tn = F.normalize(target, dim=-1, eps=1e-8)
        cos = torch.sum(pn * tn, dim=-1).clamp(-1 + 1e-7, 1 - 1e-7)
        return torch.acos(cos) * (180.0 / math.pi)  # (B,M)

    # ---------------- 训练一步 ----------------
    def train_step(self, data, target, face_idx=None):
        self.model.train()
        if isinstance(data, np.ndarray):
            data = torch.from_numpy(data)
        if isinstance(target, np.ndarray):
            target = torch.from_numpy(target)
        data = data.to(self.device, non_blocking=True)
        target = target.to(self.device, non_blocking=True)

        # 前处理：旋转(+X) + L2（并返回 R、Y_rot 以兼容旧路径）
        with torch.no_grad():
            data_rot, y_glb, R, y_rot = self._preprocess(data, target)

        self.opt.zero_grad(set_to_none=True)
        out = self.model(
            data_rot,
            patch_grad_s=self.patch_grad_scale,
            detach_patch=self.patch_detach,
        )

        
        if isinstance(face_idx, np.ndarray):
            face_idx = torch.from_numpy(face_idx)
        face_idx = face_idx.to(self.device, non_blocking=True).long()
        loss, metrics = self.compute_loss_two_head_faceagg(out, y_glb, R, face_idx)

        loss.backward()

        total_norm = torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)


        self.opt.step()

        if self.logger is not None and "loss_main" in metrics:
            self.logger.add_scalar("train/loss_main", metrics["loss_main"])
            self.logger.add_scalar("train/loss_face", metrics["loss_face"])
        return metrics["loss"]

    # ---------------- 验证（MSE，旋转系；若要也做面级聚合，可加 face_idx 并复用 faceagg） ----------------
    def evaluate(self, dev_loader, sampling_size):
        self.model.eval()
        total_loss, total_cnt = 0.0, 0
        with torch.no_grad():
            for i in range(dev_loader.length()):
                data_batches, label_batches = dev_loader.generate_batch(i, sampling_size)
                for j in range(data_batches.shape[0]):
                    d = torch.from_numpy(data_batches[j]).to(self.device)
                    y = torch.from_numpy(label_batches[j]).to(self.device)

                    # 前处理：得到旋转到局部系的数据与标签
                    d_rot, y_glb, R, y_rot = self._preprocess(d, y)

                    out = self.model(d_rot, detach_patch=True)
                    n_hat = out[0] if isinstance(out, tuple) else out
                    n_hat = F.normalize(n_hat, dim=-1, eps=1e-8)
                    y_n   = F.normalize(y_rot, dim=-1, eps=1e-8)
                    loss  = F.mse_loss(n_hat, y_n)
                    B = d.shape[0]
                    total_loss += float(loss.detach().cpu()) * B
                    total_cnt  += B
        return total_loss / max(1, total_cnt)

    # ---------------- 验证角度（度） ----------------
    def evaluate_angle_faceagg(self, dev_loader, sampling_size):
        """
        面级聚合评估（角度，单位：度）：
        1) 先把每个 patch 的预测从局部旋回到全局；
        2) 按 face_idx 聚合同一面的多次预测为一个向量（平均）；
        3) 与全局 GT 对齐后计算夹角（每个 face 一次），最后取均值。
        """
        self.model.eval()
        total_ang_sum, total_face_cnt = 0.0, 0
        with torch.no_grad():
            for i in range(dev_loader.length()):
                # 需要 dev loader 返回三元：(data, label, face_idx)
                ret = dev_loader.generate_batch(i, sampling_size)
                if len(ret) != 3:
                    raise RuntimeError("dev_loader must be constructed with return_face_idx=True for face-agg eval.")
                data_batches, label_batches, idx_batches = ret

                for j in range(data_batches.shape[0]):
                    d  = torch.from_numpy(data_batches[j]).to(self.device)
                    y  = torch.from_numpy(label_batches[j]).to(self.device)
                    ib = torch.from_numpy(idx_batches[j]).to(self.device).long()  # (B,M)

                    # 前处理：旋到局部系，但保留全局GT与旋转矩阵R
                    d_rot, y_glb, R, _y_rot = self._preprocess(d, y)
                    out = self.model(d_rot, detach_patch=True)
                    n_hat_rot = out[0] if isinstance(out, tuple) else out  # (B,M,3) in rotated coords

                    # 反旋回到全局
                    n_hat_glb = self._unrotate_vec(n_hat_rot, R)  # (B,M,3)

                    # 按面聚合：同一 face 的多次预测向量求平均
                    pred_face = self._agg_by_face(n_hat_glb, ib)  # (U,3)
                    gt_face   = self._agg_by_face(y_glb,   ib)    # (U,3)  # y_glb 本就是全局GT

                    # 计算每个 face 的角度（度）
                    pn  = torch.nn.functional.normalize(pred_face, dim=-1, eps=1e-8)
                    tn  = torch.nn.functional.normalize(gt_face,   dim=-1, eps=1e-8)
                    cos = (pn * tn).sum(dim=-1).clamp(-1 + 1e-7, 1 - 1e-7)  # (U,)
                    ang = torch.acos(cos) * (180.0 / math.pi)               # (U,)

                    # 以“face 等权”口径累计
                    U = ang.numel()
                    total_ang_sum += float(ang.sum().detach().cpu())
                    total_face_cnt += int(U)

        return total_ang_sum / max(1, total_face_cnt)
