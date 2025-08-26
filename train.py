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
from torch.optim.lr_scheduler import ReduceLROnPlateau
from tqdm import tqdm

from train_utils.fileloader import Loader
from train_utils.checkpoints import CheckpointIO
from train_utils.trainer import Trainer
from train_utils import config

# =======================================================
# 数据路径（按你的目录调整）
# =======================================================
DATASET_ROOT_TRAIN = "dataset_s_i1/train1001"
DATASET_ROOT_DEV   = "dataset_s_i1/dev1001"
PATCH_ROOT         = "patches"
META_TRAIN         = os.path.join(DATASET_ROOT_TRAIN, "meta.json")
META_DEV           = os.path.join(DATASET_ROOT_DEV,   "meta.json")

# =======================================================
# 运行开关 & 恢复设置
# =======================================================
RUN_S1 = False                   # 是否执行 S1（只训 Face-Encoder + 线性头）
RUN_S2 = True                   # 是否执行 S2（只训 Patch-Encoder）

RESUME_S1 = True                # 是否尝试从 stage1_latest.pt 恢复（含优化器/调度器）
RESUME_S2 = True                # 是否尝试从 stage2_latest.pt 恢复（含优化器/调度器）

# 当 S2 无法从 latest 恢复时，是否在开始前加载 S1 的 best 权重
LOAD_S1_BEST_BEFORE_S2 = True

# =======================================================
# Scheduler & Early Stop（两阶段可独立配置）
# =======================================================
# —— S1 ——
S1_EPOCHS      = 10
S1_BATCH_SIZE  = 8
S1_LR          = 1e-4
USE_LR_SCHED_S1        = True      # 启用 ReduceLROnPlateau
USE_EARLY_STOP_S1      = True      # 启用早停
S1_SCHED_FACTOR        = 0.5
S1_SCHED_PATIENCE      = 2
S1_MIN_LR              = 1e-6
S1_EARLY_STOP_PATIENCE = 3

# —— S2 ——
S2_EPOCHS      = 20
S2_BATCH_SIZE  = 8
S2_LR          = 1e-4
USE_LR_SCHED_S2        = True
USE_EARLY_STOP_S2      = True
S2_SCHED_FACTOR        = 0.5
S2_SCHED_PATIENCE      = 2
S2_MIN_LR              = 1e-6
S2_EARLY_STOP_PATIENCE = 5   # 放宽，避免与调度器撞车

# 早停与“有提升”的容差
IMPROVE_DELTA         = 2e-6

# =======================================================
# 训练细节
# =======================================================
SHUFFLE_TRAIN  = True
DROP_LAST      = True
LOG_CENTER_ONLY_ANGLE = False  # 如需打印 center-only 角度，设 True

# =======================================================
# 日志/输出（支持固定目录 / 环境变量覆盖 / 追加日志）
# =======================================================
# 如果你想复用同一个目录（便于 resume），把 USE_FIXED_OUT_DIR 设为 True 并设置 FIXED_OUT_DIR；
# 也可以通过环境变量 OUT_DIR 覆盖（优先级更高）。
USE_FIXED_OUT_DIR = True
FIXED_OUT_DIR     = os.path.join('out', 'two_stage_1001')  # 自行修改为你的固定实验目录
APPEND_LOGS       = True  # 复用目录时追加日志而不是覆盖

OUT_DIR_ENV = os.getenv('OUT_DIR')
if OUT_DIR_ENV and len(OUT_DIR_ENV.strip()) > 0:
    out_dir = OUT_DIR_ENV
elif USE_FIXED_OUT_DIR:
    out_dir = FIXED_OUT_DIR
else:
    current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    out_dir = os.path.join('out', current_time)

os.makedirs(out_dir, exist_ok=True)

logger = SummaryWriter(os.path.join(out_dir, 'logs'))
log_mode = 'a' if APPEND_LOGS else 'w'
logfile = open(os.path.join(out_dir, 'log.txt'), log_mode, buffering=1)

def log_print(msg: str):
    print(msg, flush=True)
    logfile.write(msg + "\n"); logfile.flush(); logger.flush()

def get_lr(optimizer):
    return optimizer.param_groups[0]['lr'] if optimizer.param_groups else float('nan')

# =======================================================
# Loader（两套，分别用不同 B）
# =======================================================
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
    drop_last=True,
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
    drop_last=True,
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

# =======================================================
# 设备 & 基础模型
# =======================================================
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
base = config.get_model(
    device,
    patch_num=patch_num,
    lsd_r_size=lsd_r_size,
    lsd_t_size=lsd_t_size,
)

# =======================================================
# 两个轻封装（共享 base）
# =======================================================
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

# =======================================================
# S1：只训 Face + 线性头
# =======================================================
def run_stage1():
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

    opt1 = optim.Adam(
        (p for p in stage1.parameters() if p.requires_grad),
        lr=S1_LR
    )
    sched1 = None
    if USE_LR_SCHED_S1:
        sched1 = ReduceLROnPlateau(opt1, mode='min', factor=S1_SCHED_FACTOR,
                                   patience=S1_SCHED_PATIENCE, min_lr=S1_MIN_LR)

    trainer1 = Trainer(stage1, opt1, device=device)
    ckpt1 = CheckpointIO(out_dir, model=stage1, optimizer=opt1, scheduler=sched1)

    # ===== 恢复（latest） =====
    start_epoch = 0
    if RESUME_S1:
        try:
            scalars = ckpt1.load('stage1_latest.pt')
            start_epoch = int(scalars.get('epoch_it', -1)) + 1
            log_print(f"[S1] Resume from epoch {start_epoch}")
        except Exception as e:
            log_print(f"[S1] Resume skipped ({e})")

    best_val = float('inf')
    no_improve = 0

    log_print("==== Stage 1: Train Face-Encoder + Linear Head (rotated patches) ====")
    for epoch in range(start_epoch, S1_EPOCHS):
        # —— 训练 ——
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

        # —— 验证 ——
        val_loss = trainer1.evaluate(dev_loader_s1, sampling_size)
        val_ang  = trainer1.evaluate_angle(dev_loader_s1, sampling_size, center_only=False)

        # —— 调度（先 step，再观察 LR 变化）——
        prev_lr = get_lr(opt1)
        if sched1 is not None:
            sched1.step(val_loss)
        new_lr = get_lr(opt1)
        if new_lr < prev_lr:
            log_print(f"[S1] LR reduced: {prev_lr:.6g} -> {new_lr:.6g} (patience hit)")

        # —— 日志 ——
        log_line = f"[S1] lr={get_lr(opt1):.6g} | val_mse: {val_loss:.6f} | angle_all: {val_ang:.3f}°"
        if LOG_CENTER_ONLY_ANGLE:
            val_ctr = trainer1.evaluate_angle(dev_loader_s1, sampling_size, center_only=True)
            log_line += f" | angle_center: {val_ctr:.3f}°"
        log_print(log_line)

        logger.add_scalar('stage1/val_loss', val_loss, epoch)
        logger.add_scalar('stage1/val_angle_deg', val_ang, epoch)

        # —— 保存 latest/best + 早停 ——
        ckpt1.save('stage1_latest.pt', epoch_it=epoch, metric=val_loss)
        if val_loss < best_val - IMPROVE_DELTA:
            best_val = val_loss; no_improve = 0
            ckpt1.save('stage1_best.pt', epoch_it=epoch, metric=best_val)
            log_print(f"[S1] Now best val loss: {best_val:.6f}")
        else:
            # 若刚刚降了学习率，给小学习率一个机会：本轮不计入无提升
            if new_lr < prev_lr:
                no_improve = 0
            else:
                no_improve += 1
            if USE_EARLY_STOP_S1 and no_improve >= S1_EARLY_STOP_PATIENCE:
                log_print("[S1] Early stop triggered.")
                break

    return stage1, best_val

# =======================================================
# S2：只训 Patch（Face 冻结）
# =======================================================
def run_stage2(resume_first=True):
    # 冻结 Face、启用 Patch
    for name, p in base.named_parameters():
        if name.startswith("face_encoder") or "face" in name:
            p.requires_grad = False
        else:
            p.requires_grad = True

    stage2 = Stage2PatchOnly(base, freeze_face=True).to(device)
    if torch.cuda.device_count() >= 2:
        stage2 = nn.DataParallel(stage2, device_ids=list(range(torch.cuda.device_count())))
        log_print(f"[S2] Using GPUs: {stage2.device_ids}")

    opt2 = optim.Adam(
        (p for p in stage2.parameters() if p.requires_grad),
        lr=S2_LR
    )
    sched2 = None
    if USE_LR_SCHED_S2:
        sched2 = ReduceLROnPlateau(opt2, mode='min', factor=S2_SCHED_FACTOR,
                                   patience=S2_SCHED_PATIENCE, min_lr=S2_MIN_LR)

    trainer2 = Trainer(stage2, opt2, device=device)
    ckpt2 = CheckpointIO(out_dir, model=stage2, optimizer=opt2, scheduler=sched2)

    # ===== 优先尝试恢复（latest） =====
    start_epoch = 0
    resumed = False
    if RESUME_S2 and resume_first:
        try:
            scalars = ckpt2.load('stage2_latest.pt')
            start_epoch = int(scalars.get('epoch_it', -1)) + 1
            resumed = True
            log_print(f"[S2] Resume from epoch {start_epoch}")
        except Exception as e:
            log_print(f"[S2] Resume skipped ({e})")

    # ===== 若未恢复，则尽量加载 S1 的 best 权重（作为 S2 起点） =====
    if (not resumed) and LOAD_S1_BEST_BEFORE_S2:
        try:
            _tmp_ckpt = CheckpointIO(out_dir, model=stage2)
            _tmp_ckpt.load('stage1_best.pt')
            log_print('[S2] Loaded stage1_best.pt into model before starting S2.')
        except Exception as e:
            log_print(f'[S2] stage1_best.pt not found or failed to load ({e}); proceed with current weights.')

    best_val = float('inf')
    no_improve = 0

    log_print("==== Stage 2: Freeze Face, Train Patch-Encoder ====")
    for epoch in range(start_epoch, S2_EPOCHS):
        # —— 训练 ——
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

        # —— 验证 ——
        val_loss = trainer2.evaluate(dev_loader_s2, sampling_size)
        val_ang  = trainer2.evaluate_angle(dev_loader_s2, sampling_size, center_only=False)

        # —— 调度（先 step，再观察 LR 变化）——
        prev_lr = get_lr(opt2)
        if sched2 is not None:
            sched2.step(val_loss)
        new_lr = get_lr(opt2)
        if new_lr < prev_lr:
            log_print(f"[S2] LR reduced: {prev_lr:.6g} -> {new_lr:.6g} (patience hit)")

        # —— 日志 ——
        log_line = f"[S2] lr={get_lr(opt2):.6g} | val_mse: {val_loss:.6f} | angle_all: {val_ang:.3f}°"
        if LOG_CENTER_ONLY_ANGLE:
            val_ctr = trainer2.evaluate_angle(dev_loader_s2, sampling_size, center_only=True)
            log_line += f" | angle_center: {val_ctr:.3f}°"
        log_print(log_line)

        logger.add_scalar('stage2/val_loss', val_loss, epoch)
        logger.add_scalar('stage2/val_angle_deg', val_ang, epoch)

        # —— 保存 latest/best + 早停 ——
        ckpt2.save('stage2_latest.pt', epoch_it=epoch, metric=val_loss)
        if val_loss < best_val - IMPROVE_DELTA:
            best_val = val_loss; no_improve = 0
            ckpt2.save('stage2_best.pt', epoch_it=epoch, metric=best_val)
            log_print(f"[S2] Now best val loss: {best_val:.6f}")
        else:
            # 若刚刚降了学习率，给小学习率一个机会：本轮不计入无提升
            if new_lr < prev_lr:
                no_improve = 0
            else:
                no_improve += 1
            if USE_EARLY_STOP_S2 and no_improve >= S2_EARLY_STOP_PATIENCE:
                log_print("[S2] Early stop triggered.")
                break

    return stage2, best_val

# =======================================================
# 主流程
# =======================================================
def main():
    if RUN_S1:
        run_stage1()
    else:
        log_print("[Main] Skip Stage 1.")

    if RUN_S2:
        run_stage2(resume_first=True)
    else:
        log_print("[Main] Skip Stage 2.")

    logger.close()
    logfile.close()

if __name__ == "__main__":
    main()
