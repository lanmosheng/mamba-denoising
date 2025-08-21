import torch
import torch.nn.functional as F
import numpy as np
import os

class Trainer:
    '''Trainer object for the Mamba-based model.
    
    Args:
        model (nn.Module): The model to be trained.
        optimizer (torch.optim.Optimizer): The optimizer to be used.
        device (torch.device): The device (GPU/CPU) for computation.
        input_type (str): The type of input data (optional, for flexibility).
        vis_dir (str): Directory for storing visualizations (optional).
        threshold (float): Threshold for classification tasks (not used here).
        eval_sample (bool): Flag to evaluate on samples during training (optional).
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
        '''Performs a training step.'''
        self.model.train()  # Set model to training mode
        self.optimizer.zero_grad()  # Reset the gradients
        
        # Move data to the device (GPU or CPU)
        data = data.to(self.device)
        label = label.to(self.device)

        # Compute the loss for the current batch
        loss = self.compute_loss(data, label)
        
        # Backpropagation
        loss.backward()
        
        # Update the model weights
        self.optimizer.step()
        
        return loss.item()  # Return the loss value for logging or monitoring

    def eval_step(self, data, label):
        '''Performs an evaluation step (without backpropagation).'''
        self.model.eval()  # Set model to evaluation mode
        
        # Move data to the device
        data = data.to(self.device)
        label = label.to(self.device)

        with torch.no_grad():  # Disable gradient computation for evaluation
            loss = self.compute_loss(data, label)
        
        return loss  # Return the loss value for logging or monitoring

    def compute_loss(self, data, label, save_debug=True, debug_prefix="debug"):
        '''Compute the loss for a batch.'''
        device = self.device

        # Move data to device and ensure it's in float32 type
        data = data if isinstance(data, torch.Tensor) else torch.as_tensor(data)
        label = label if isinstance(label, torch.Tensor) else torch.as_tensor(label)
        
        data = data.to(device=device, dtype=torch.float32)
        label = label.to(device=device, dtype=torch.float32)

        # Forward pass: Get the model's output
        output = self.model.pred(data)  # Expected shape: [B, M, 3] (predicted normals)

        # Basic sanity check for NaN or Inf values
        def _check(name, t):
            if torch.isnan(t).any() or torch.isinf(t).any():
                if save_debug:
                    torch.save(
                        {"output": output.detach().cpu(),
                         "data": data.detach().cpu(),
                         "label": label.detach().cpu()},
                        f"{debug_prefix}_dump.pt"
                    )
                raise ValueError(f"[Loss Debug] {name} contains NaN/Inf")

        _check("model output", output)
        _check("label", label)

        # Core: Use normalized MSE (mean squared error) loss on unit vectors
        output_n = F.normalize(output, dim=1, eps=1e-8)
        label_n = F.normalize(label, dim=1, eps=1e-8)
        loss = F.mse_loss(output_n, label_n)  # MSE loss between predicted and true normalized normals
        
        if torch.isnan(loss) or torch.isinf(loss):
            if save_debug:
                torch.save(
                    {"output_n": output_n.detach().cpu(),
                     "label_n": label_n.detach().cpu()},
                    f"{debug_prefix}_norm_dump.pt"
                )
            raise ValueError("[Loss Debug] Loss is NaN/Inf")

        return loss  # Return the computed loss
