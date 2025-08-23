# -*- coding: utf-8 -*-
import os
from tqdm import tqdm
# —— 指定两张卡（需在 import torch 之前）——
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"
# 可选：减少碎片
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "max_split_size_mb:128")

import time
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter

from train_utils.fileloader import Loader
from train_utils.checkpoints import CheckpointIO
from train_utils.trainer import Trainer
from train_utils import config

# 允许 TF32（提速，显存影响小）
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

# ----------------- 基本配置 -----------------
batch_size = 10
epochs = 20
current_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
out_dir = os.path.join('out', current_time)
os.makedirs(out_dir, exist_ok=True)

# ----------------- 数据加载器 -----------------
train_loader = Loader(
    dataset_root="dataset_s_i1/train1001/",
    patch_root="patches/",
    batch_size=batch_size,
    meta_path="dataset_s_i1/train1001/meta.json",
    drop_last=True,
    shuffle_faces=True,
    mmap=True,
    rotation_anchor='center',  # or 'mean'
)

dev_loader = Loader(
    dataset_root="dataset_s_i1/dev1001/",
    patch_root="patches/",
    batch_size=batch_size,
    meta_path="dataset_s_i1/dev1001/meta.json",
    drop_last=False,
    shuffle_faces=False,
    mmap=True,
    rotation_anchor='center',
)

meta = train_loader.get_meta_data()
patch_num  = meta['patch_num']
lsd_r_size = meta['lsd_r_size']
lsd_t_size = meta['lsd_t_size']
sampling_size = lsd_r_size * lsd_t_size + 1

# ----------------- 设备 & 模型 -----------------
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
model = config.get_model(
    device,
    patch_num=patch_num,
    lsd_r_size=lsd_r_size,
    lsd_t_size=lsd_t_size,
    # 可按需改：d_model、层数、位置编码等
)

if torch.cuda.device_count() >= 2:
    model = nn.DataParallel(model, device_ids=[0, 1])
    print(f"Using GPUs: {model.device_ids}", flush=True)
else:
    print("Single GPU or CPU.", flush=True)

optimizer = optim.Adam(model.parameters(), lr=1e-4)

# ----------------- Trainer / 日志 / 检查点 -----------------
trainer = Trainer(
    model, optimizer, device=device,
    use_amp=True, amp_dtype=torch.float16,
    microbatch_size=None,   # 如显存仍紧张可设 8/16 做梯度累积
)

logger = SummaryWriter(os.path.join(out_dir, 'logs'))
logfile = open(os.path.join(out_dir, 'log.txt'), 'w', buffering=1)  # 行缓冲

checkpoint_io = CheckpointIO(out_dir, model=model, optimizer=optimizer)

# ----------------- 训练主循环 -----------------
metric_val_best = float('inf')

for epoch_it in range(epochs):
    for i in range(train_loader.length()):
        total_batches = train_loader.count_batches(i)  # 需要在 Loader 里实现这个函数
        pbar = tqdm(total=total_batches, desc=f"Epoch {epoch_it} | mesh {i}", leave=False)

        mesh_loss_sum = 0.0
        mesh_sample_cnt = 0

        for tdata, tlabel in train_loader.iter_batches(i, sampling_size):
            cur_B = tdata.shape[0]
            loss = trainer.train_step(tdata, tlabel)
            logger.add_scalar('train/loss_batch', loss, epoch_it); logger.flush()

            mesh_loss_sum  += float(loss) * cur_B
            mesh_sample_cnt += cur_B
            pbar.update(1)

        pbar.close()

        if mesh_sample_cnt > 0:
            mesh_avg_loss = mesh_loss_sum / mesh_sample_cnt
            print(f"[Epoch {epoch_it:02d}] file={i:03d}: avg_loss={mesh_avg_loss:.6f}", flush=True)
            logfile.write(f"[Epoch {epoch_it:02d}] file={i:03d}: avg_loss={mesh_avg_loss:.6f}\n"); logfile.flush()
            logger.add_scalar('train/loss_mesh', mesh_avg_loss, epoch_it); logger.flush()


    # 保存最新检查点
    checkpoint_io.save('model.pt', epoch_it=epoch_it)

    # 验证
    metric_val = trainer.evaluate(dev_loader, sampling_size)
    logfile.write(f"Validation metric : {metric_val:.6f}\n")
    logfile.flush()

    if metric_val < metric_val_best:
        metric_val_best = metric_val
        checkpoint_io.save('model_best.pt', epoch_it=epoch_it)
    logfile.write(f"Now best model metric: {metric_val_best:.6f}\n")

logger.close()
logfile.close()
