# -*- coding: utf-8 -*-
import torch

try:
    from train_utils.model_mamba import TwoStageMamba
except ImportError:
    # 若你放在 coder.py，则改为 from coder import TwoStageMamba
    from train_utils.model_mamba import TwoStageMamba

def get_model(
    device: torch.device,
    *,
    patch_num: int,
    lsd_r_size: int,
    lsd_t_size: int,
    d_model: int = 64,
    face_depth: int = 4,
    patch_depth: int = 4,
    add_pos_emb: bool = True,     # LSD token（N维序列）的极坐标位置编码
    add_patch_pos: bool = False,  # patch 内位置编码（可选）
    d_state: int = 16,
    d_conv: int = 4,
    expand: int = 2,
    residual_scale: float = 0.5,
    p_drop: float = 0.05,
):
    """
    构建两级 Mamba 模型：
      - FaceEncoder:  (B, M, N, 3) -> (B, M, d)
      - PatchEncoder: (B, M, d)    -> (B, M, 3)
    其中 N = 1 + lsd_r_size * lsd_t_size，M = patch_num。
    """
    model = TwoStageMamba(
        patch_num=patch_num,
        lsd_r_size=lsd_r_size,
        lsd_t_size=lsd_t_size,
        d_model=d_model,
        face_depth=face_depth,
        patch_depth=patch_depth,
        add_pos_emb=add_pos_emb,
        add_patch_pos=add_patch_pos,
        d_state=d_state,
        d_conv=d_conv,
        expand=expand,
        residual_scale=residual_scale,
        p_drop=p_drop,
    )
    return model.to(device)
