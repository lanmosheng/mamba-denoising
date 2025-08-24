# -*- coding: utf-8 -*-
import os
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0,1")
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "max_split_size_mb:128")

import time
import contextlib
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm

from train_utils.fileloader import Loader
from train_utils.checkpoints import CheckpointIO
from train_utils.trainer import Trainer
from train_utils import config

# =========================
# 路径配置（按你的目录改）
# =========================
DATASET_ROOT_TRAIN = "dataset_s_i1/train1001"
DATASET_ROOT_DEV   = "dataset_s_i1/dev1001"
PATCH_ROOT         = "patches"
META_TRAIN         = os.path.join(DATASET_ROOT_TRAIN, "meta.json")
META_DEV           = os.path.join(DATASET_ROOT_DEV,   "meta.json")

# =========================
# 阶段 1（只训 Face-Encoder + 线性头）
# =========================
S1_EPOCHS      = 10
S1_BATCH_SIZE  = 8
S1_LR          = 1e-4

# =========================
# 阶段 2（只训 Patch-Encoder）
# =========================
S2_EPOCHS      = 10
S2_BATCH_SIZE  = 8
S2_LR          = 1e-4

SHUFFLE_TRAIN  = True
DROP_LAST      = True

# =========================
# 日志/输出
# =========================
current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
out_dir = os.path.join('out', current_time)
os.makedirs(out_dir, exist_ok=True)
logger = SummaryWriter(os.path.join(out_dir, 'logs'))
logfile = open(os.path.join(out_dir, 'log.txt'), 'w', buffering=1)

def log_print(msg: str):
    print(msg, flush=True)
    logfile.write(msg + "\n"); logfile.flush(); logger.flush()

# =========================
# Loader（两套，分别用不同 B）
# 说明：Loader 内已做“统一旋转到 +X + patch级 z-score”，S1/S2 输入分布一致
# =========================
train_loader_s1 = Loader(
    dataset_root=DATASET_ROOT_TRAIN,
    patch_root=PATCH_ROOT,
    batch_size=S1_BATCH_SIZE,
    meta_path=META_TRAIN,
    drop_last=DROP_LAST,
    shuffle_faces=SHUFFLE_TRAIN,
    mmap=True,
    rotation_anchor='center',
)

dev_loader_s1 = Loader(
    dataset_root=DATASET_ROOT_DEV,
    patch_root=PATCH_ROOT,
    batch_size=S1_BATCH_SIZE,
    meta_path=META_DEV,
    drop_last=False,
    shuffle_faces=False,
    mmap=True,
    rotation_anchor='center',
)

train_loader_s2 = Loader(
    dataset_root=DATASET_ROOT_TRAIN,
    patch_root=PATCH_ROOT,
    batch_size=S2_BATCH_SIZE,
    meta_path=META_TRAIN,
    drop_last=DROP_LAST,
    shuffle_faces=SHUFFLE_TRAIN,
    mmap=True,
    rotation_anchor='center',
)

dev_loader_s2 = Loader(
    dataset_root=DATASET_ROOT_DEV,
    patch_root=PATCH_ROOT,
    batch_size=S2_BATCH_SIZE,
    meta_path=META_DEV,
    drop_last=False,
    shuffle_faces=False,
    mmap=True,
    rotation_anchor='center',
)

# 形状信息
meta = train_loader_s1.get_meta_data()
patch_num  = meta['patch_num']
lsd_r_size = meta['lsd_r_size']
lsd_t_size = meta['lsd_t_size']
sampling_size = lsd_r_size * lsd_t_size + 1

# =========================
# 设备 & 基础模型
# =========================
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
base = config.get_model(
    device,
    patch_num=patch_num,
    lsd_r_size=lsd_r_size,
    lsd_t_size=lsd_t_size,
)

# =========================
# 两个轻封装（共享 base）
# forward 接受 **kwargs（兼容 Trainer 在评估时传入的 detach_patch 等）
# =========================
class Stage1FaceOnly(nn.Module):
    """S1：Face-Encoder + 线性头 (d->3)，输出 (B,M,3)。"""
    def __init__(self, base_model: nn.Module):
        super().__init__()
        self.base = base_model
        d_guess = getattr(self.base, 'd_model', None)
        if d_guess is None:
            d_guess = getattr(getattr(self.base, 'face', object()), 'd_model', 64)
        self.face_head = nn.Linear(int(d_guess), 3)
    def forward(self, X, **kwargs):
        F = self.base.face_encoder(X)              # (B,M,d)
        y = self.face_head(F)                      # (B,M,3)
        return torch.nn.functional.normalize(y, dim=-1, eps=1e-8)

class Stage2PatchOnly(nn.Module):
    """S2：冻结 Face（no_grad 取 F），只训 Patch-Encoder，输出 (B,M,3)。"""
    def __init__(self, base_model: nn.Module, freeze_face: bool = True):
        super().__init__()
        self.base = base_model
        self.freeze_face = freeze_face
    def forward(self, X, **kwargs):
        ctx = torch.no_grad() if self.freeze_face else contextlib.nullcontext()
        with ctx:
            F = self.base.face_encoder(X)          # (B,M,d)
        y = self.base.patch_encoder(F)             # (B,M,3)
        return y

# =========================
# S1：只训 Face + 线性头
# =========================
# 冻结 Patch、启用 Face
for name, p in base.named_parameters():
    if name.startswith("patch_encoder") or "patch" in name:
        p.requires_grad = False
    else:
        p.requires_grad = True

stage1 = Stage1FaceOnly(base).to(device)
if torch.cuda.device_count() >= 2:
    stage1 = nn.DataParallel(stage1, device_ids=list(range(torch.cuda.device_count())))
    log_print(f"[S1] Using GPUs: {stage1.device_ids}")

opt1 = optim.Adam((p for p in stage1.parameters() if p.requires_grad), lr=S1_LR)
trainer1 = Trainer(stage1, opt1, device=device)  # ← 与当前 trainer.py 构造签名匹配
ckpt1 = CheckpointIO(out_dir, model=stage1, optimizer=opt1)

best_val_1 = float('inf')
log_print("==== Stage 1: Train Face-Encoder + Linear Head (rotated patches) ====")
for epoch in range(S1_EPOCHS):
    # 训练
    for mi in range(train_loader_s1.length()):
        total_batches = train_loader_s1.count_batches(mi)
        pbar = tqdm(total=total_batches, desc=f"[S1] Epoch {epoch} | mesh {mi}", leave=False)
        mesh_loss_sum, mesh_cnt = 0.0, 0
        for Xb, Yb in train_loader_s1.iter_batches(mi, sampling_size):
            B_cur = Xb.shape[0]
            loss = trainer1.train_step(Xb, Yb)
            mesh_loss_sum += float(loss) * B_cur
            mesh_cnt      += B_cur
            pbar.update(1)
        pbar.close()
        if mesh_cnt:
            mesh_avg = mesh_loss_sum / mesh_cnt
            log_print(f"[S1][Epoch {epoch:02d}] file={mi:03d}: avg_loss={mesh_avg:.6f}")
            logger.add_scalar('stage1/train_loss_mesh', mesh_avg, epoch)
    # 验证（loss + angle）
    val_loss = trainer1.evaluate(dev_loader_s1, sampling_size)
    val_ang  = trainer1.evaluate_angle(dev_loader_s1, sampling_size, center_only=False)
    log_print(f"[S1] Validation loss: {val_loss:.6f} | angle: {val_ang:.3f}°")
    logger.add_scalar('stage1/val_loss', val_loss, epoch)
    logger.add_scalar('stage1/val_angle_deg', val_ang, epoch)
    ckpt1.save('stage1_latest.pt', epoch_it=epoch, metric=val_loss)
    if val_loss < best_val_1:
        best_val_1 = val_loss
        ckpt1.save('stage1_best.pt', epoch_it=epoch, metric=best_val_1)
        log_print(f"[S1] Now best val loss: {best_val_1:.6f}")

# =========================
# S2：只训 Patch（Face 冻结）
# =========================
for name, p in base.named_parameters():
    if name.startswith("face_encoder") or "face" in name:
        p.requires_grad = False
    else:
        p.requires_grad = True

stage2 = Stage2PatchOnly(base, freeze_face=True).to(device)
if torch.cuda.device_count() >= 2:
    stage2 = nn.DataParallel(stage2, device_ids=list(range(torch.cuda.device_count())))
    log_print(f"[S2] Using GPUs: {stage2.device_ids}")

opt2 = optim.Adam((p for p in stage2.parameters() if p.requires_grad), lr=S2_LR)
trainer2 = Trainer(stage2, opt2, device=device)  # ← 与当前 trainer.py 构造签名匹配
ckpt2 = CheckpointIO(out_dir, model=stage2, optimizer=opt2)

best_val_2 = float('inf')
log_print("==== Stage 2: Freeze Face, Train Patch-Encoder ====")
for epoch in range(S2_EPOCHS):
    # 训练
    for mi in range(train_loader_s2.length()):
        total_batches = train_loader_s2.count_batches(mi)
        pbar = tqdm(total=total_batches, desc=f"[S2] Epoch {epoch} | mesh {mi}", leave=False)
        mesh_loss_sum, mesh_cnt = 0.0, 0
        for Xb, Yb in train_loader_s2.iter_batches(mi, sampling_size):
            B_cur = Xb.shape[0]
            loss = trainer2.train_step(Xb, Yb)
            mesh_loss_sum += float(loss) * B_cur
            mesh_cnt      += B_cur
            pbar.update(1)
        pbar.close()
        if mesh_cnt:
            mesh_avg = mesh_loss_sum / mesh_cnt
            log_print(f"[S2][Epoch {epoch:02d}] file={mi:03d}: avg_loss={mesh_avg:.6f}")
            logger.add_scalar('stage2/train_loss_mesh', mesh_avg, epoch)
    # 验证（loss + angle）
    val_loss = trainer2.evaluate(dev_loader_s2, sampling_size)
    val_ang  = trainer2.evaluate_angle(dev_loader_s2, sampling_size, center_only=False)
    log_print(f"[S2] Validation loss: {val_loss:.6f} | angle: {val_ang:.3f}°")
    logger.add_scalar('stage2/val_loss', val_loss, epoch)
    logger.add_scalar('stage2/val_angle_deg', val_ang, epoch)
    ckpt2.save('stage2_latest.pt', epoch_it=epoch, metric=val_loss)
    if val_loss < best_val_2:
        best_val_2 = val_loss
        ckpt2.save('stage2_best.pt', epoch_it=epoch, metric=best_val_2)
        log_print(f"[S2] Now best val loss: {best_val_2:.6f}")

logger.close()
logfile.close()
