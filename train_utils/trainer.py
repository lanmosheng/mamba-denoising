import os
from tqdm import tqdm
from tqdm import trange
import torch
from torch.nn import functional as F
from torch import distributions as dist
import numpy


class Trainer():
    ''' Trainer object for the Occupancy Network.

    Args:
        model (nn.Module): Occupancy Network model
        optimizer (optimizer): pytorch optimizer object
        device (device): pytorch device
        input_type (str): input type
        vis_dir (str): visualization directory
        threshold (float): threshold value
        eval_sample (bool): whether to evaluate samples

    '''



    def __init__(self, model, optimizer, device=None, input_type='img',
                 vis_dir=None, threshold=0.5, eval_sample=False):
        self.model = model
        self.optimizer = optimizer
        self.device = device
        self.input_type = input_type
        self.vis_dir = vis_dir
        self.threshold = threshold
        self.eval_sample = eval_sample

        if vis_dir is not None and not os.path.exists(vis_dir):
            os.makedirs(vis_dir)

    def train_step(self, data, label):
        ''' Performs a training step.

        Args:
            data (dict): data dictionary
        '''
        self.model.train()
        self.optimizer.zero_grad()

        loss= self.compute_loss(data, label)
        loss.backward()
        self.optimizer.step()
        return loss.item()
    
    
    def evaluate(self, val_loader, sampling_size):
        ''' Performs an evaluation.
        Args:
            val_loader (dataloader): pytorch dataloader
        '''
        val=0.0
        num=0
        for i in range(val_loader.length()):
            tdata, tlabel=val_loader.generate_batch(i, sampling_size)
            for j in range(tdata.shape[0]):
                loss = self.eval_step(tdata[j], tlabel[j])
                val=val+torch.sum(loss.float())
                num=num+tdata[j].shape[0]

        return val/num
    
    def eval_step(self, data, label):
        ''' Performs an evaluation step.

        Args:
            data (dict): data dictionary
        '''
        self.model.eval()
        
        data =torch.tensor(data).to(self.device).float()
        label = torch.tensor(label).to(self.device).float()        

        with torch.no_grad():
            loss = self.model.compute_loss(data, label)
        return  loss

    # def compute_loss(self, data, label):
    #     ''' Computes the loss.

    #     Args:
    #         data (dict): data dictionary
    #     '''
    #     device = self.device
    #     data =torch.tensor(data).to(self.device).float()
    #     label = torch.tensor(label).to(self.device).float()   
    #     output = self.model.pred(data)
    #     loss_fn = torch.nn.MSELoss()
    #     loss = loss_fn(output, label)

    #     return loss.float()

    # def compute_loss(self, data, label):
    #     device = self.device
    #     data = torch.tensor(data).to(device).float()
    #     label = torch.tensor(label).to(device).float()

    #     output = self.model.pred(data)

    #     # ✅ 检查 output 和 label 是否有异常
    #     if torch.isnan(output).any() or torch.isinf(output).any():
    #         print("❌ [Loss Debug] Model output contains NaN or Inf")
    #         torch.save(output, "debug_output.pt")
    #         torch.save(data, "debug_input.pt")
    #         torch.save(label, "debug_label.pt")
    #         raise ValueError("Model output contains NaN or Inf")

    #     loss_fn = torch.nn.MSELoss()
    #     loss = loss_fn(output, label)

    #     if torch.isnan(loss).any() or torch.isinf(loss).any():
    #         print("❌ [Loss Debug] Loss value is NaN or Inf")
    #         raise ValueError("Loss value is NaN or Inf")

    #     return loss.float()
    def compute_loss(self, data, label, *, save_debug=True, debug_prefix="debug"):
        device = self.device

        # 更稳的张量转换（若已是 Tensor 不额外复制）
        data  = data  if isinstance(data,  torch.Tensor) else torch.as_tensor(data)
        label = label if isinstance(label, torch.Tensor) else torch.as_tensor(label)
        data  = data.to(device=device, dtype=torch.float32)
        label = label.to(device=device, dtype=torch.float32)

        # 前向
        output = self.model.pred(data)  # 期望形状: [B, 3]

        # 基本健诊
        def _check(name, t):
            if torch.isnan(t).any() or torch.isinf(t).any():
                if save_debug:
                    torch.save(
                        {"output": output.detach().cpu(),
                        "data":   data.detach().cpu(),
                        "label":  label.detach().cpu()},
                        f"{debug_prefix}_dump.pt"
                    )
                raise ValueError(f"[Loss Debug] {name} contains NaN/Inf")

        _check("model output", output)
        _check("label", label)
        # 如需也检查输入，取消下一行注释
        # _check("input data", data)

        # —— 核心：单位化后的 MSE（与 1-cos 等价，更稳）——
        output_n = F.normalize(output, dim=1, eps=1e-8)
        label_n  = F.normalize(label,  dim=1, eps=1e-8)
        loss = F.mse_loss(output_n, label_n)

        if torch.isnan(loss) or torch.isinf(loss):
            if save_debug:
                torch.save(
                    {"output_n": output_n.detach().cpu(),
                    "label_n":  label_n.detach().cpu()},
                    f"{debug_prefix}_norm_dump.pt"
                )
            raise ValueError("[Loss Debug] Loss is NaN/Inf")

        return loss


    