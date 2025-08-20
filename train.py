import sys
import torch
import torch.optim as optim
import numpy as np
import os
os.environ["CUDA_VISIBLE_DEVICES"] = '0'  # 单卡：只用 0 号 GPU
from tensorboardX import SummaryWriter
import argparse
import time
import config
from trainer import Trainer
from checkpoints import CheckpointIO
from fileloader import Loader
import pickle
import random
from mamba_ssm import Mamba
import torch.nn as nn
from torch.optim.lr_scheduler import ReduceLROnPlateau
import math
from torch.optim.lr_scheduler import LambdaLR

#---------------------------------------------------#
#   设置种子
#---------------------------------------------------#
def seed_everything(seed=11):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

if __name__ == '__main__':  
    lsd_r = int(sys.argv[1])
    lsd_t = int(sys.argv[2])
    sampling_size = lsd_r * lsd_t + 1
    is_cuda = (torch.cuda.is_available() )
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)
    seed_everything(0)
    t0 = time.time()

    # ---- 模型 → to(device) ----
    model = config.get_model(device, sampling_size=sampling_size, lsd_t=lsd_t, lsd_r=lsd_r)
    model = model.to(device)

    # ---- Optimizer / Trainer ----
    base_lr = 1e-4     # 目标最大学习率（warmup 结束时达到）
    min_lr  = 3e-6     # 余弦衰减下限（可改 1e-6 ~ 3e-6）
    warmup_epochs = 2  # warmup 轮数
    total_epochs = 20
    optimizer = optim.Adam(model.parameters(), lr=base_lr)
    # ---- 目录先创建，再打开日志，再建 writer ----
    # out_dir = f'out{sampling_size}/'
    # out_dir = f'out_s2/'
    out_dir = 'out_1601_drop'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    run_tag = time.strftime('%Y%m%d-%H%M%S')
    logfile = open(os.path.join(out_dir, f'log_{run_tag}.txt'), 'w', buffering=1)
    logger = SummaryWriter(os.path.join(out_dir, 'logs'))

    # ---- CheckpointIO 绑定当前模型 ----
    checkpoint_io = CheckpointIO(out_dir, model=model, optimizer=optimizer)
    try:
        load_dict = checkpoint_io.load('model.pt')
        print("Resumed training from checkpoint.")
    except FileNotFoundError:
        load_dict = dict()
        print("No checkpoint found. Starting training from scratch.")
    epoch_it = load_dict.get('epoch_it', -1)
    it = load_dict.get('it', -1)   # 训练步计数器
    metric_val_best = load_dict.get('loss_val_best', np.inf) 
    

    # ---- 定义 Warmup + Cosine（在 load 之后构建；对齐 last_epoch）----
    def lr_lambda(epoch):
        if epoch < warmup_epochs:
            # 线性：min_lr -> base_lr
            scale = (epoch + 1) / max(1, warmup_epochs)
            target = min_lr + (base_lr - min_lr) * scale
        else:
            # 余弦阶段同理，加 +1 对齐下一轮
            prog = (epoch - warmup_epochs + 1) / max(1, total_epochs - warmup_epochs)
            prog = min(max(prog, 0.0), 1.0)  # 夹一下，防越界
            target = min_lr + 0.5 * (base_lr - min_lr) * (1.0 + math.cos(math.pi * prog))
        return target / base_lr

    scheduler = LambdaLR(optimizer, lr_lambda=lr_lambda, last_epoch=epoch_it)

    if epoch_it == -1:
        for pg in optimizer.param_groups:
            pg['lr'] = min_lr

    trainer = Trainer(model, optimizer, device=device)


    batch_size = 160

    nparameters = sum(p.numel() for p in model.parameters())
    logfile.write('Total number of parameters: %d\n' % nparameters)

    print("pos emb on? ->", getattr(model.predictor, "pos_emb", None) is not None)
    #test_loader=Loader("test/",batch_size)
    dev_loader   = Loader(f'dev{sampling_size}/',   batch_size)
    train_loader = Loader(f'train{sampling_size}/', batch_size)


    for epoch_it in range(epoch_it + 1, total_epochs):
        logfile.flush()
        cur_lr = optimizer.param_groups[0]['lr']
        logger.add_scalar('lr', float(cur_lr), epoch_it)
        logfile.write(f'[Epoch {epoch_it:02d}] lr={cur_lr:.6e}\n')
        for i in range(train_loader.length()):
            tdata, tlabel = train_loader.generate_batch(i, sampling_size=sampling_size)
            # print(tdata.shape)
            file_loss_sum = 0.0
            for j in range(tdata.shape[0]):
                loss = trainer.train_step(tdata[j], tlabel[j])
                it += 1
                logger.add_scalar('train/loss', float(loss), it)
                file_loss_sum += float(loss)
            file_loss = file_loss_sum / max(1, tdata.shape[0])
            logger.add_scalar('train/loss_file_avg', file_loss, it)
            print('[Epoch %02d] file=%02d: loss=%.6f' % (epoch_it, i, file_loss))
            logfile.write('[Epoch %02d] file=%02d: loss=%.6f\n' % (epoch_it, i, file_loss))

        logfile.write('Saving checkpoint\n')
        metric_val2 = trainer.evaluate(dev_loader, sampling_size=sampling_size)
        metric_val2 = float(metric_val2)  # 转为 float 以便格式化
        logfile.write('Validation metric : %.6f\n' % metric_val2)
        
        
        scheduler.step()
        cur_lr_after = optimizer.param_groups[0]['lr']
        logfile.write('Learning rate : %.6e\n' % cur_lr_after)

        checkpoint_io.save('model.pt', epoch_it=epoch_it, it=it, loss_val_best=metric_val_best)
        if metric_val2 < metric_val_best:
            metric_val_best = metric_val2
            logfile.write('New best model (loss %.6f)\n' % metric_val_best)
            checkpoint_io.save('model_best.pt', epoch_it=epoch_it, it=it, loss_val_best=metric_val_best)

        checkpoint_io.save('model_'+str(epoch_it)+'.pt', epoch_it=epoch_it, it=it, loss_val_best=metric_val_best)
    logfile.write('Training Over, best model (loss %.6f)\n' % metric_val_best)
    logger.close()
    logfile.close()
