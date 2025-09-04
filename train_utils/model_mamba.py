# -*- coding: utf-8 -*-
import math
from typing import Optional
import torch
import torch.nn as nn
import torch.nn.functional as F
from mamba_ssm import Mamba
from typing import Tuple

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
    输入:  X ∈ [B, N, 3]  (单个面的 LSD，N=1001 等)
    输出:  n ∈ [B, 3]     (单位法向量)
    """
    def __init__(self,
                 N: int,
                 d_model: int = 64,
                 depth: int = 4,
                 d_state: int = 16,
                 d_conv: int = 4,
                 expand: int = 2,
                 residual_scale: float = 0.5,
                 p_drop: float = 0.05):
        super().__init__()
        self.N = N
        self.d_model = d_model

        # 仅保留基础投影
        self.input_proj = nn.Linear(3, d_model, bias=True)

        # Mamba 堆叠
        self.layers = nn.ModuleList([
            MambaBlock(d_model, d_state=d_state, d_conv=d_conv, expand=expand,
                       residual_scale=residual_scale, p_drop=p_drop)
            for _ in range(depth)
        ])

        # 池化 + 归一化 + 法向预测头
        self.out_norm = RMSNorm(d_model)
        self.head = nn.Linear(d_model, 3, bias=True)


    def forward(self, X: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        X: [B, N, 3]
        return:
        n    ∈ [B, 3]        # 单位法向量
        feat ∈ [B, d_model]  # 面级 embedding（供 S2 / 离线缓存）
        """
        B, N, C = X.shape
        assert C == 3, f"FaceEncoder: 最后维应为3，实得 {C}"
        assert N == self.N, f"FaceEncoder: 期望 N={self.N}, 实得 {N}"

        h = self.input_proj(X)    # [B, N, d]
        for blk in self.layers:
            h = blk(h)            # [B, N, d]

        h    = h.mean(dim=1)      # [B, d]
        feat = self.out_norm(h)   # [B, d]
        n    = self.head(feat)    # [B, 3]
        n    = F.normalize(n, dim=1)
        return n, feat

    def pred(self, X: torch.Tensor) -> torch.Tensor:
        """X: [B, N, 3] -> n: [B, 3]"""
        self.eval()
        with torch.no_grad():
            n, _ = self.forward(X)
        return n
    
    def extract_feat(self, X: torch.Tensor) -> torch.Tensor:
        """X: [B, N, 3] -> feat: [B, d_model]"""
        self.eval()
        with torch.no_grad():
            _, feat = self.forward(X)
        return feat


# ---------------- 第二层：Patch-Encoder (Mamba₂) ----------------
class PatchEncoder(nn.Module):
    """
    输入:  F ∈ [B, M, d_model]
    输出:  Y ∈ [B, M, 3]（单位化法向量）
    """
    def __init__(self,
                 M: int,
                 d_model: int = 64,
                 depth: int = 4,
                 d_state: int = 16,
                 d_conv: int = 4,
                 expand: int = 2,
                 residual_scale: float = 0.5,
                 p_drop: float = 0.05):
        super().__init__()
        self.M = M
        self.d_model = d_model

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
        assert d == self.d_model, f"PatchEncoder: 期望特征维 d_model={self.d_model}, 实得 {d}"

        h = Fm
        if self.pos_emb_patch is not None:
            idx = torch.arange(M, device=Fm.device).unsqueeze(0).expand(B, M)
            h = h + self.pos_emb_patch(idx)

        for blk in self.layers:
            h = blk(h)

        h = self.head_norm(h)
        y = self.head(h)              # (B, M, 3)
        return F.normalize(y, dim=-1, eps=1e-8)

    # PatchEncoder: forward 已返回单位化法向 [B, M, 3]；pred 直接复用
    def pred(self, Fm: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        """Fm: [B, M, d_model] -> y: [B, M, 3]（可选 mask: [B, M]）"""
        self.eval()
        with torch.no_grad():
            y = self.forward(Fm)  # 已单位化
            if mask is not None:
                # 将无效位（-1处）置零，便于后续聚合时跳过
                y = y.masked_fill(~mask[..., None], 0.0)
        return y
