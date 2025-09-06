# train_s1_mamba.py
# -*- coding: utf-8 -*-
import os, time, math, random
from typing import Optional
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1")
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "max_split_size_mb:128")
from tqdm import tqdm
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.tensorboard import SummaryWriter

# =========================
# 路径 & 实验输出
# =========================
DATASET_ROOT_TRAIN = "dataset_s_i1/train1601"
DATASET_ROOT_DEV   = "dataset_s_i1/dev1601"
PATCH_ROOT         = "train_patches"  # 若 S1FaceLoaderV2 内部使用 _resolve_patch_path，会用到
META_TRAIN         = os.path.join(DATASET_ROOT_TRAIN, "meta.json")
META_DEV           = os.path.join(DATASET_ROOT_DEV,   "meta.json")

OUT_NAME           = "s1_mamba_face_1601"
OUT_DIR            = os.path.join("out", OUT_NAME)
APPEND_LOGS        = True  # 追加日志而不是覆盖

# =========================
# 训练超参（S1）
# =========================
EPOCHS            = 100
BATCH_SIZE_TRAIN  = 240
BATCH_SIZE_DEV    = 240
LR                = 1e-4
WEIGHT_DECAY      = 0.0
GRAD_CLIP_NORM    = 1.0   # None 关闭

# 学习率调度
USE_LR_SCHED      = True
SCHED_FACTOR      = 0.5
SCHED_PATIENCE    = 3
MIN_LR            = 1e-6

# =========================
# Loader 设置（S1 专用）
# =========================
SLICE_RATIO_TRAIN = 0.10  # 每个 epoch 每个 mesh 取 10% patch，跨 epoch 不重叠
SLICE_RATIO_DEV   = 1.00  # 验证用全量
SHUFFLE_TRAIN     = True
DROP_LAST_TRAIN   = True
SHUFFLE_DEV       = False
DROP_LAST_DEV     = False
SHUFFLE_SEED      = 0     # 复现性

# =========================
# 运行设备 & 随机数
# =========================
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
SEED   = 42

# =========================
# 导入你项目内的模块（全部使用新版）
# =========================
from train_utils.fileloader import S1FaceLoaderV2   # 仅旋转规范化：中心面 token0 → (1,0,0)
from train_utils import config                           # 内含 get_face_mamba(...)

# ============== 工具函数 ==============
def seed_everything(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def angle_deg(pred: torch.Tensor, target: torch.Tensor, eps: float = 1e-8) -> torch.Tensor:
    """pred/target: (...,3) -> 每样本角度(度)"""
    pn = F.normalize(pred, dim=-1, eps=eps)
    tn = F.normalize(target, dim=-1, eps=eps)
    cos = torch.sum(pn * tn, dim=-1).clamp(-1.0 + 1e-7, 1.0 - 1e-7)
    return torch.acos(cos) * (180.0 / math.pi)

def get_out_dir() -> str:
    os.makedirs(OUT_DIR, exist_ok=True)
    return OUT_DIR

# ============== 训练/验证 ==============
def train_one_epoch(model, optimizer, loader, epoch, logger, grad_clip_norm):
    model.train()
    loader.set_epoch(epoch)
    total_batches = loader.total_batches()
    refresh_every = max(1, total_batches // 100)
    total_loss, total_cnt = 0.0, 0
    # 用 tqdm 包装 batch 生成器；total 不设也可以，tqdm 会动态推进
    pbar = tqdm(
        loader.iter_batches(),
        total=(total_batches if total_batches > 0 else None),
        desc=f"Epoch {epoch}",
        unit="batch",
        dynamic_ncols=True,
        mininterval=0.5,
        smoothing=0.1,
    )

    for step, batch in enumerate(pbar, 1):
        X = torch.from_numpy(batch['X']).to(DEVICE, non_blocking=True)
        Y = torch.from_numpy(batch['Y']).to(DEVICE, non_blocking=True)

        optimizer.zero_grad(set_to_none=True)
        out = model(X)
        n_pred = out[0] if isinstance(out, tuple) else out

        n_pred_n = F.normalize(n_pred, dim=-1)
        Y_n      = F.normalize(Y,      dim=-1)
        loss = F.mse_loss(n_pred_n, Y_n)

        loss.backward()
        if grad_clip_norm is not None:
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
        optimizer.step()

        B = X.shape[0]
        total_loss += float(loss.detach().cpu()) * B
        total_cnt  += B

        # 进度条上显示当前 loss（同时也可显示 lr）

        if logger is not None:
            logger.add_scalar("train_s1/mse", float(loss.detach().cpu()))

    pbar.close()  
    return total_loss / max(1, total_cnt)

@torch.no_grad()
def validate_angle(model: nn.Module, loader: S1FaceLoaderV2, logger: Optional[SummaryWriter]) -> float:
    model.eval()
    ang_sum, ang_cnt = 0.0, 0
    for batch in loader.iter_batches():
        X = torch.from_numpy(batch['X']).to(DEVICE, non_blocking=True)
        Y = torch.from_numpy(batch['Y']).to(DEVICE, non_blocking=True)
        out = model(X)
        n_pred = out[0] if isinstance(out, tuple) else out
        ang = angle_deg(n_pred, Y)   # (B,)
        ang_sum += float(ang.sum().cpu())
        ang_cnt += int(ang.numel())
    mean_deg = ang_sum / max(1, ang_cnt)
    if logger is not None:
        logger.add_scalar("val_s1/angle_deg", mean_deg)
    return mean_deg

# ============== 主流程 ==============
def main():
    seed_everything(SEED)
    out_dir = get_out_dir()
    log_mode = 'a' if APPEND_LOGS else 'w'
    logfile = open(os.path.join(out_dir, "log.txt"), log_mode, buffering=1, encoding="utf-8")
    def log_print(msg: str):
        print(msg, flush=True)
        logfile.write(msg + "\n"); logfile.flush()

    logger = SummaryWriter(os.path.join(out_dir, "logs"))

    # ---------- Loader ----------
    train_loader = S1FaceLoaderV2(
        dataset_root=DATASET_ROOT_TRAIN,
        batch_size=BATCH_SIZE_TRAIN,
        meta_path=META_TRAIN,
        patch_root=PATCH_ROOT,
        slice_ratio=SLICE_RATIO_TRAIN,
        shuffle_faces=SHUFFLE_TRAIN,
        shuffle_seed=SHUFFLE_SEED,
        drop_last=DROP_LAST_TRAIN,
        mmap_mode='r',
    )
    dev_loader = S1FaceLoaderV2(
        dataset_root=DATASET_ROOT_DEV,
        batch_size=BATCH_SIZE_DEV,
        meta_path=META_DEV,
        patch_root=PATCH_ROOT,
        slice_ratio=SLICE_RATIO_DEV,     # 全量评估
        shuffle_faces=SHUFFLE_DEV,
        shuffle_seed=SHUFFLE_SEED,
        drop_last=DROP_LAST_DEV,
        mmap_mode='r',
    )

    # 读取 meta（用于构建模型 & 打印）
    meta = train_loader.meta
    lsd_r_size = int(meta['lsd_r_size'])
    lsd_t_size = int(meta['lsd_t_size'])
    sampling_size = 1 + lsd_r_size * lsd_t_size
    log_print(f"[Meta] r={lsd_r_size}, t={lsd_t_size}, N={sampling_size}")

    # ---------- Model（仅新版） ----------
    face_model = config.get_face_mamba(
        DEVICE,
        lsd_r_size=lsd_r_size,
        lsd_t_size=lsd_t_size,
        d_model=64,
        depth=4,
        d_state=16,
        d_conv=4,
        expand=2,
        residual_scale=0.5,
        p_drop=0.05,
    )
    log_print("[Model] FaceEncoder via config.get_face_mamba(...)")

    # 多卡
    if torch.cuda.device_count() >= 2:
        face_model = nn.DataParallel(face_model, device_ids=list(range(torch.cuda.device_count())))
        log_print(f"[S1] Using GPUs: {face_model.device_ids}")

    # ---------- Optim & Sched ----------
    opt = optim.AdamW(
        (p for p in face_model.parameters() if p.requires_grad),
        lr=LR, weight_decay=WEIGHT_DECAY
    )
    sched = ReduceLROnPlateau(opt, mode='min', factor=SCHED_FACTOR,
                              patience=SCHED_PATIENCE, min_lr=MIN_LR) if USE_LR_SCHED else None

    # ---------- Resume ----------
    latest_path = os.path.join(out_dir, "stage1_latest.pt")
    best_path   = os.path.join(out_dir, "stage1_best.pt")
    start_epoch = 0
    best_val = float('inf')
    if os.path.exists(best_path):
        try:
            obj = torch.load(best_path, map_location='cpu')
            (face_model.module if isinstance(face_model, nn.DataParallel) else face_model)\
                .load_state_dict(obj['model'], strict=False)
            opt.load_state_dict(obj['optimizer'])
            start_epoch = int(obj.get('epoch', -1)) + 1
            best_val = float(obj.get('best_val', float('inf')))
            log_print(f"[Resume] stage1_best.pt from epoch {start_epoch}")
        except Exception as e:
            log_print(f"[Resume] skip ({e})")

    # ---------- Train Loop ----------
    log_print("==== Train Stage 1: FaceEncoder (MSE train / angle val) ====")
    for epoch in range(start_epoch, EPOCHS):
        train_mse = train_one_epoch(face_model, opt, train_loader, epoch, logger, GRAD_CLIP_NORM)
        val_deg   = validate_angle(face_model, dev_loader, logger)

        prev_lr = opt.param_groups[0]['lr']
        if sched is not None:
            sched.step(val_deg)
        cur_lr = opt.param_groups[0]['lr']
        if cur_lr < prev_lr:
            log_print(f"[S1] LR reduced: {prev_lr:.3e} -> {cur_lr:.3e}")

        log_print(f"[S1][Epoch {epoch:03d}] train_mse={train_mse:.6f} | val_angle={val_deg:.3f}° | lr={cur_lr:.2e}")

        # 保存 latest
        state = {
            "model": (face_model.module if isinstance(face_model, nn.DataParallel) else face_model).state_dict(),
            "optimizer": opt.state_dict(),
            "epoch": epoch,
            "val_angle": val_deg,
            "best_val": best_val,
        }
        torch.save(state, latest_path)

        # 保存 best
        if val_deg < best_val - 1e-6:
            best_val = val_deg
            torch.save(state, best_path)
            log_print(f"[S1] New best val angle: {best_val:.6f}°")

    logger.close()
    logfile.close()

if __name__ == "__main__":
    main()
