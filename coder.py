import torch.nn as nn
import torch
import torch.nn.functional as F
import numpy as np
from mamba_ssm import Mamba

class Network(nn.Module):

    def __init__(self, predictor, device):
        super().__init__()

        self.predictor = predictor.to(device)
        self._device = device



    def compute_loss(self, data, lab):
        output = F.normalize(self.pred(data), dim=1)
        loss = torch.nn.functional.cosine_similarity(output, lab, dim=1)
        loss2=torch.acos(loss)*180.0/np.pi
        isnan=torch.isnan(loss2)
        for i in range(isnan.shape[0]):
            if isnan[i]==True:
                loss2[i]=0
        return loss2

    def pred(self, data):
        output = self.predictor(data)
        return output


# ---- 基础：更稳的归一化（建议替代 BatchNorm）----
class RMSNorm(nn.Module):
    def __init__(self, d_model: int, eps: float = 1e-8):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(d_model))
        self.eps = eps
    def forward(self, x):
        # x: [..., d_model]
        norm = x.pow(2).mean(dim=-1, keepdim=True).add(self.eps).rsqrt()
        return self.weight * x * norm

# ---- Mamba 块：Pre-Norm + 残差缩放（更稳、不易 NaN）----
class MambaBlock(nn.Module):
    def __init__(self, d_model: int, d_state: int = 16, d_conv: int = 4,
                 expand: int = 2, residual_scale: float = 0.5, p_drop: float = 0.0):
        super().__init__()
        self.norm = RMSNorm(d_model)  # Pre-Norm（图里的 N）
        self.core = Mamba(d_model=d_model, d_state=d_state, d_conv=d_conv, expand=expand)
        self.drop = nn.Dropout(p_drop)
        self.res_scale = residual_scale

    def forward(self, x):             # x: [B, L, C]
        h = self.core(self.norm(x))   # 先归一化再进 Mamba
        h = self.drop(h)
        return x + self.res_scale * h # 残差 + 缩放

# ---- 最终模型 ----
class MambaNet(nn.Module):
    def __init__(self,
                 num_patches: int = 1601,          # = 1 + lsd_r_size * lsd_t_size
                 lsd_r_size: int | None = None,
                 lsd_t_size: int | None = None,
                 d_model: int = 96,
                 d_state: int = 16,
                 d_conv: int = 4,
                 expand: int = 2,
                 depth: int = 8,
                 residual_scale: float = 0.5,
                 p_drop: float = 0.0,
                 add_pos_emb: bool = True):
        super().__init__()
        self.num_patches = num_patches
        self.add_pos_emb = add_pos_emb

        self.input_proj = nn.Linear(3, d_model)

        if add_pos_emb:
            assert lsd_r_size is not None and lsd_t_size is not None, \
                "Please pass lsd_r_size and lsd_t_size when add_pos_emb=True"
            assert num_patches == 1 + lsd_r_size * lsd_t_size, \
                f"L must be 1 + r*t, got L={num_patches}, r={lsd_r_size}, t={lsd_t_size}"
            self.lsd_r_size = lsd_r_size
            self.lsd_t_size = lsd_t_size

            self.pos_emb_ring = nn.Embedding(lsd_r_size + 1, d_model)
            self.pos_emb_ang  = nn.Embedding(lsd_t_size,     d_model)
            self.polar_proj   = nn.Linear(3, d_model, bias=False)

            import math, torch as _torch
            L = num_patches
            t = _torch.arange(L, dtype=_torch.long)
            ring_idx = _torch.zeros(L, dtype=_torch.long)
            ang_idx  = _torch.zeros(L, dtype=_torch.long)
            theta    = _torch.zeros(L, dtype=_torch.float)
            r_norm   = _torch.zeros(L, dtype=_torch.float)
            ring_idx[0] = 0; ang_idx[0] = 0; theta[0] = 0.0; r_norm[0] = 0.0
            if L > 1:
                t1 = t[1:]
                ring = 1 + (t1 - 1) // lsd_t_size
                ang  =      (t1 - 1) %  lsd_t_size
                ring_idx[1:] = ring
                ang_idx[1:]  = ang
                theta[1:]    = 2.0 * math.pi * ang.float() / float(lsd_t_size)
                r_norm[1:]   = ring.float() / float(lsd_r_size)
            self.register_buffer("ring_idx_buf", ring_idx, persistent=False)
            self.register_buffer("ang_idx_buf",  ang_idx,  persistent=False)
            polar_feat = _torch.stack([_torch.cos(theta), _torch.sin(theta), r_norm], dim=-1)  # [L,3]
            self.register_buffer("polar_feat_buf", polar_feat, persistent=False)
        else:
            self.register_buffer("ring_idx_buf", torch.empty(0, dtype=torch.long), persistent=False)
            self.register_buffer("ang_idx_buf",  torch.empty(0, dtype=torch.long), persistent=False)
            self.register_buffer("polar_feat_buf", torch.empty(0, 3), persistent=False)

        self.layers = nn.ModuleList([
            MambaBlock(d_model, d_state, d_conv, expand,
                      residual_scale=residual_scale, p_drop=p_drop)
            for _ in range(depth)
        ])

        self.head_norm = RMSNorm(d_model)
        self.output_head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.GELU(),
            nn.Linear(128, 3),
        )

    def forward(self, x):                        # x: [B, L, 3]
        B, L, Cin = x.shape
        assert L == self.num_patches, f"Expected L={self.num_patches}, got {L}"

        x = self.input_proj(x)                   # [B, L, d_model]

        if self.add_pos_emb:
            pe = self.pos_emb_ring(self.ring_idx_buf) \
               + self.pos_emb_ang(self.ang_idx_buf) \
               + self.polar_proj(self.polar_feat_buf)     # [L, d_model]
            x = x + pe.unsqueeze(0)                        # [B, L, d_model]

        for blk in self.layers:
            x = blk(x)                           # [B, L, d_model]

        x = x.mean(dim=1)                        # [B, d_model]
        x = self.head_norm(x)
        y = self.output_head(x)                  # [B, 3]
        y = F.normalize(y, dim=1, eps=1e-8)
        return y

# class MambaNet(nn.Module):
#     def __init__(self,
#                  num_patches: int = 1001,
#                  d_model: int = 96,
#                  d_state: int = 16,
#                  d_conv: int = 4,
#                  expand: int = 2,
#                  depth: int = 6,
#                  residual_scale: float = 0.5,
#                  p_drop: float = 0.0,
#                  add_pos_emb: bool = True):
#         """
#         输入: x ∈ [B, num_patches, 3]
#         输出: y ∈ [B, 3]（单位向量）
#         """
#         super().__init__()
#         self.num_patches = num_patches
#         self.add_pos_emb = add_pos_emb

#         self.input_proj = nn.Linear(3, d_model)

#         if add_pos_emb:
#             # 可学习位置编码（当 1001 个采样点有固定顺序时建议打开）
#             self.pos_emb = nn.Embedding(num_patches, d_model)
#         else:
#             self.pos_emb = None

#         self.layers = nn.ModuleList([
#             MambaBlock(d_model, d_state, d_conv, expand,
#                       residual_scale=residual_scale, p_drop=p_drop)
#             for _ in range(depth)
#         ])

#         # Head：无 BatchNorm/Tanh，先做归一化再 MLP，最后单位化输出
#         self.head_norm = RMSNorm(d_model)
#         self.output_head = nn.Sequential(
#             nn.Linear(d_model, 128),
#             nn.GELU(),
#             nn.Linear(128, 3),
#         )

#     def forward(self, x):                        # x: [B, L, 3]
#         B, L, Cin = x.shape
#         assert L == self.num_patches, f"Expected L={self.num_patches}, got {L}"

#         x = self.input_proj(x)                   # [B, L, d_model]

#         if self.pos_emb is not None:
#             pos = torch.arange(L, device=x.device).unsqueeze(0).expand(B, L)  # [B, L]
#             x = x + self.pos_emb(pos)           # 加位置

#         for blk in self.layers:
#             x = blk(x)                           # Mamba 层堆叠（Pre-Norm + 残差）

#         x = x.mean(dim=1)                        # mean pooling: [B, d_model]
#         x = self.head_norm(x)                    # 稳定一点
#         y = self.output_head(x)                  # [B, 3]
#         y = F.normalize(y, dim=1, eps=1e-8)      # 单位化，替代 Tanh
#         return y




