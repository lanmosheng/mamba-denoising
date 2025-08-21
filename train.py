import torch
from tqdm import tqdm
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
from tensorboardX import SummaryWriter
import time
from train_utils.trainer import Trainer
from train_utils.checkpoints import CheckpointIO
from train_utils.fileloader import Loader
from train_utils import config

# Configurations can be moved to a separate config.py
batch_size = 80
current_time = time.strftime("%Y-%m-%d %H", time.localtime())
out_dir = f'out{current_time}/'

# Dataset
train_loader = Loader("train/", batch_size)
dev_loader = Loader("dev/", batch_size)

train_loader = Loader(
    dataset_root="dataset/train/",  # lsd.npy 和 gt.npy 所在目录
    patch_root="patches/",           # patches 目录
    batch_size=batch_size,
    meta_path="dataset/meta.json",   # 或者指定 meta 文件路径
    drop_last=True,
    shuffle_faces=True,
    rotation_anchor='center'
)

# Meta

meta = Loader.get_meta_data()
lsd_r_size, lsd_t_size, patch_num = meta['lsd_r_size'], meta['lsd_t_size'], meta['patch_num']
sampling_size = lsd_r_size * lsd_t_size + 1


# Training
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = config.get_model(device)
optimizer = optim.Adam(model.parameters(), lr=1e-4)
trainer = Trainer(model, optimizer, device=device)
# Logging
logger = SummaryWriter(os.path.join(out_dir, 'logs'))
logfile = open(os.path.join(out_dir, 'log.txt'), 'w')
checkpoint_io = CheckpointIO(out_dir, model=model, optimizer=optimizer)


# Main Training Loop
for epoch_it in range(20):
    for i in tqdm(range(train_loader.length())):
        tdata, tlabel = train_loader.generate_batch(i)
        loss = trainer.train_step(tdata, tlabel)
        logger.add_scalar('train/loss', loss, epoch_it)

        print(f"[Epoch {epoch_it}] file={i}: loss={loss:.6f}")
        logfile.write(f"[Epoch {epoch_it}] file={i}: loss={loss:.6f}\n")

    # Save model checkpoints
    checkpoint_io.save('model.pt', epoch_it=epoch_it)
    
    # Validation
    metric_val = trainer.evaluate(dev_loader)
    logfile.write(f"Validation metric : {metric_val:.6f}\n")
    
    if metric_val < metric_val_best:
        metric_val_best = metric_val
        checkpoint_io.save('model_best.pt', epoch_it=epoch_it)

# Closing logs
logger.close()
