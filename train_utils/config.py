import torch
from train_utils.model_mamba import FaceEncoder
from train_utils.model_mamba import PatchEncoder
def get_face_mamba(
    device: torch.device,
    *,
    lsd_r_size: int,
    lsd_t_size: int,
    d_model: int = 64,
    depth: int = 4,
    d_state: int = 16,
    d_conv: int = 4,
    expand: int = 2,
    residual_scale: float = 0.5,
    p_drop: float = 0.05,
):
    """
    FaceEncoder: X ∈ [B, N, 3] -> (n ∈ [B,3], feat ∈ [B, d_model])
    其中 N = 1 + lsd_r_size * lsd_t_size
    """
    N = 1 + int(lsd_r_size) * int(lsd_t_size)
    model = FaceEncoder(
        N=N,
        d_model=d_model,
        depth=depth,
        d_state=d_state,
        d_conv=d_conv,
        expand=expand,
        residual_scale=residual_scale,
        p_drop=p_drop,
    )
    return model.to(device)


def get_patch_mamba(
    device: torch.device,
    *,
    patch_num: int,
    d_model: int = 64,
    depth: int = 4,
    d_state: int = 16,
    d_conv: int = 4,
    expand: int = 2,
    residual_scale: float = 0.5,
    p_drop: float = 0.05,
):
    """
    PatchEncoder: F ∈ [B, M, d_model] -> Y ∈ [B, M, 3]
    其中 M = patch_num
    """
    model = PatchEncoder(
        M=int(patch_num),
        d_model=d_model,
        depth=depth,
        d_state=d_state,
        d_conv=d_conv,
        expand=expand,
        residual_scale=residual_scale,
        p_drop=p_drop,
    )
    return model.to(device)


def get_models(
    device: torch.device,
    *,
    patch_num: int,
    lsd_r_size: int,
    lsd_t_size: int,
    d_model: int = 64,
    face_depth: int = 4,
    patch_depth: int = 4,
    d_state: int = 16,
    d_conv: int = 4,
    expand: int = 2,
    residual_scale: float = 0.5,
    p_drop: float = 0.05,
):
    """
    便捷函数：一次性返回 (face_model, patch_model)
    """
    face_model = get_face_mamba(
        device,
        lsd_r_size=lsd_r_size,
        lsd_t_size=lsd_t_size,
        d_model=d_model,
        depth=face_depth,
        d_state=d_state,
        d_conv=d_conv,
        expand=expand,
        residual_scale=residual_scale,
        p_drop=p_drop,
    )
    patch_model = get_patch_mamba(
        device,
        patch_num=patch_num,
        d_model=d_model,
        depth=patch_depth,
        d_state=d_state,
        d_conv=d_conv,
        expand=expand,
        residual_scale=residual_scale,
        p_drop=p_drop,
    )
    return face_model, patch_model
