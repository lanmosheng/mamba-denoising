# -*- coding: utf-8 -*-
import math
from typing import Optional
import torch
import torch.nn as nn
import torch.nn.functional as F
from mamba_ssm import Mamba

# ---------------- RMSNorm / MambaBlock ----------------
class RMSNorm(nn.Module):
    def __init__(self, d_model: int, eps: float = 1e-8):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(d_model))
        self.eps = eps
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        norm = x.pow(2).mean(dim=-1, keepdim=True).add(self.eps).rsqrt()
        return self.weight * x * norm

class MambaBlock(nn.Module):
    def __init__(self, d_model: int, d_state: int = 16, d_conv: int = 4,
                 expand: int = 2, residual_scale: float = 0.5, p_drop: float = 0.0):
        super().__init__()
        self.norm = RMSNorm(d_model)
        self.core = Mamba(d_model=d_model, d_state=d_state, d_conv=d_conv, expand=expand)
        self.drop = nn.Dropout(p_drop)
        self.res_scale = residual_scale
    def forward(self, x: torch.Tensor) -> torch.Tensor:  # x: [B, L, C]
        h = self.core(self.norm(x))
        h = self.drop(h)
        return x + self.res_scale * h

# ---------------- 第一层：Face-Encoder (Mamba₁) ----------------
class FaceEncoder(nn.Module):
    """
    输入:  X ∈ [B, M, N, 3]  (B: batch, M: patch_num, N: 1 + r*t)
    输出:  F ∈ [B, M, d]
    """
    def __init__(self,
                 N: int,
                 d_model: int = 64,
                 depth: int = 4,
                 lsd_r_size: Optional[int] = None,
                 lsd_t_size: Optional[int] = None,
                 add_pos_emb: bool = True,
                 d_state: int = 16,
                 d_conv: int = 4,
                 expand: int = 2,
                 residual_scale: float = 0.5,
                 p_drop: float = 0.05):
        super().__init__()
        self.N = N
        self.d_model = d_model
        self.add_pos_emb = add_pos_emb

        self.input_proj = nn.Linear(3, d_model, bias=True)

        if add_pos_emb:
            assert lsd_r_size is not None and lsd_t_size is not None, \
                "FaceEncoder(add_pos_emb=True) 需要 lsd_r_size / lsd_t_size"
            assert N == 1 + lsd_r_size * lsd_t_size, \
                f"N 应为 1+r*t，实际 N={N}, r={lsd_r_size}, t={lsd_t_size}"

            self.pos_emb_ring = nn.Embedding(lsd_r_size + 1, d_model)
            self.pos_emb_ang  = nn.Embedding(lsd_t_size, d_model)
            self.polar_proj   = nn.Linear(3, d_model, bias=False)

            ring_idx = torch.zeros(N, dtype=torch.long)
            ang_idx  = torch.zeros(N, dtype=torch.long)
            theta    = torch.zeros(N, dtype=torch.float)
            r_norm   = torch.zeros(N, dtype=torch.float)
            ring_idx[0] = 0; ang_idx[0] = 0; theta[0] = 0.0; r_norm[0] = 0.0
            if N > 1:
                t1 = torch.arange(1, N, dtype=torch.long)
                ring = 1 + (t1 - 1) // lsd_t_size
                ang  =      (t1 - 1) %  lsd_t_size
                ring_idx[1:] = ring
                ang_idx[1:]  = ang
                theta[1:]    = 2.0 * math.pi * ang.float() / float(lsd_t_size)
                r_norm[1:]   = ring.float() / float(lsd_r_size)
            polar_feat = torch.stack([torch.cos(theta), torch.sin(theta), r_norm], dim=-1)  # [N,3]

            self.register_buffer("ring_idx_buf", ring_idx, persistent=False)
            self.register_buffer("ang_idx_buf",  ang_idx,  persistent=False)
            self.register_buffer("polar_feat_buf", polar_feat, persistent=False)
        else:
            self.register_buffer("ring_idx_buf", torch.empty(0, dtype=torch.long), persistent=False)
            self.register_buffer("ang_idx_buf",  torch.empty(0, dtype=torch.long), persistent=False)
            self.register_buffer("polar_feat_buf", torch.empty(0, 3), persistent=False)

        self.layers = nn.ModuleList([
            MambaBlock(d_model, d_state=d_state, d_conv=d_conv, expand=expand,
                       residual_scale=residual_scale, p_drop=p_drop)
            for _ in range(depth)
        ])
        self.out_norm = RMSNorm(d_model)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        B, M, N, _ = X.shape
        assert N == self.N, f"FaceEncoder: 期望 N={self.N}, 实得 {N}"

        x = X.reshape(B * M, N, 3)          # (B*M, N, 3)
        h = self.input_proj(x)              # (B*M, N, d)

        if self.add_pos_emb:
            pe = self.pos_emb_ring(self.ring_idx_buf) \
               + self.pos_emb_ang(self.ang_idx_buf) \
               + self.polar_proj(self.polar_feat_buf)  # [N, d]
            h = h + pe.unsqueeze(0)

        for blk in self.layers:
            h = blk(h)                      # (B*M, N, d)

        h = h.mean(dim=1)                   # (B*M, d)
        h = self.out_norm(h)
        return h.reshape(B, M, self.d_model)

# ---------------- 第二层：Patch-Encoder (Mamba₂) ----------------
class PatchEncoder(nn.Module):
    """
    输入:  F ∈ [B, M, d]
    输出:  Y ∈ [B, M, 3]（单位化法向量）
    """
    def __init__(self,
                 M: int,
                 d_model: int = 64,
                 depth: int = 4,
                 add_patch_pos: bool = False,
                 d_state: int = 16,
                 d_conv: int = 4,
                 expand: int = 2,
                 residual_scale: float = 0.5,
                 p_drop: float = 0.05,
                 in_dim: Optional[int] = None):
        super().__init__()
        self.M = M
        self.d_model = d_model
        self.pos_emb_patch = nn.Embedding(M, d_model) if add_patch_pos else None

        self.in_dim = in_dim if in_dim is not None else d_model
        self.in_proj = nn.Linear(self.in_dim, d_model) if self.in_dim != d_model else nn.Identity()

        self.layers = nn.ModuleList([
            MambaBlock(d_model, d_state=d_state, d_conv=d_conv, expand=expand,
                       residual_scale=residual_scale, p_drop=p_drop)
            for _ in range(depth)
        ])
        self.head_norm = RMSNorm(d_model)
        self.head = nn.Linear(d_model, 3)

    def forward(self, Fm: torch.Tensor) -> torch.Tensor:
        B, M, d = Fm.shape
        assert M == self.M, f"PatchEncoder: 期望 M={self.M}, 实得 {M}"
        h = self.in_proj(Fm)
        if self.pos_emb_patch is not None:
            idx = torch.arange(M, device=Fm.device).unsqueeze(0).expand(B, M)
            h = h + self.pos_emb_patch(idx)
        for blk in self.layers:
            h = blk(h)
        h = self.head_norm(h)
        y = self.head(h)                    # (B, M, 3)
        return F.normalize(y, dim=-1, eps=1e-8)

# ---------------- 两级封装：TwoStageMamba ----------------
def grad_scale(x: torch.Tensor, s: float) -> torch.Tensor:
    """Return x with only s-fraction of gradient flowing back to its producers."""
    if s <= 0.0:
        return x.detach()
    if s >= 1.0:
        return x
    return (x - x.detach()) * s + x.detach()


class TwoStageMamba(nn.Module):
    def __init__(
        self,
        patch_num: int,              # M
        lsd_r_size: int,             # r
        lsd_t_size: int,             # t
        d_model: int = 64,
        face_depth: int = 4,
        patch_depth: int = 4,
        add_pos_emb: bool = True,    # 面内极坐标位置编码
        add_patch_pos: bool = False, # patch 内位置编码
        d_state: int = 16,
        d_conv: int = 4,
        expand: int = 2,
        residual_scale: float = 0.5,
        p_drop: float = 0.05,
        # 以下是你实验特性（可保留）
        residual_final: bool = True,
        grad_scale_s: float = 0.0,
        small_dim: int = 32,
    ):
        super().__init__()
        # 暴露给外部（S1 线性头读取）
        self.d_model = d_model
        self.M = int(patch_num)
        self.N = 1 + int(lsd_r_size) * int(lsd_t_size)

        self.face_encoder = FaceEncoder(
            N=self.N,
            d_model=d_model,
            depth=face_depth,
            lsd_r_size=lsd_r_size,
            lsd_t_size=lsd_t_size,
            add_pos_emb=add_pos_emb,
            d_state=d_state,
            d_conv=d_conv,
            expand=expand,
            residual_scale=residual_scale,
            p_drop=p_drop,
        )
        self.face_small_proj = nn.Linear(d_model, small_dim)

        in_dim = 3 + small_dim
        self.patch_encoder = PatchEncoder(
            M=self.M,
            d_model=d_model,
            depth=patch_depth,
            add_patch_pos=add_patch_pos,
            d_state=d_state,
            d_conv=d_conv,
            expand=expand,
            residual_scale=residual_scale,
            p_drop=p_drop,
            in_dim=in_dim
        )

        # ↓ 这些是你实验开关；若暂不用，也可以先关掉避免干扰
        self.residual_final = residual_final
        self.grad_scale_s = float(grad_scale_s)

        # 新：一层的“法向初稿”头（d -> 3）

        self.face_aux_head = nn.Linear(d_model, 3)

        # 你原来可能还有的模块/参数初始化...

    def _face_encode_chunked(self, X: torch.Tensor) -> torch.Tensor:
        # ← 保持你原来的实现（按 M 分块等）
        return self.face_encoder(X)

    def forward(
        self,
        X: torch.Tensor,
        patch_grad_s: float | None = None,   # 新：临时覆盖梯度比例
        detach_patch: bool = False,          # 新：强制完全切断（等价 patch_grad_s=0）
    ):
        """
        X: (B, M, N, 3)
        返回：
          - (n_hat, n1)    # 主输出 + 一层初稿（用于辅助损失）
        """
        F_face = self._face_encode_chunked(X)    # (B, M, d)

        # 一层初稿 n1（用于辅助监督 & 残差式最终输出）
        n1 = F.normalize(self.face_aux_head(F_face), dim=-1, eps=1e-8)  # (B,M,3)


        # 控制二层梯度回流到一层的比例
        s = 0.0 if detach_patch else (self.grad_scale_s if patch_grad_s is None else float(patch_grad_s))

        # 二层输出（默认当作残差 Δ）

        f_small = self.face_small_proj(F_face)
        z_in = torch.cat([n1, f_small], dim=-1)
        z_in = grad_scale(z_in, s)
        patch_out = self.patch_encoder(z_in)

        if self.residual_final:
            n_hat = F.normalize(n1 + patch_out, dim=-1, eps=1e-8)
        else:
            n_hat = F.normalize(patch_out, dim=-1, eps=1e-8)

        return (n_hat, n1)

    # 可选：推理时只要主输出
    def pred(self, X: torch.Tensor) -> torch.Tensor:
        out = self.forward(X)
        return out[0] if isinstance(out, tuple) else out