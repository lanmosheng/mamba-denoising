import torch

def stat_tensor(tensor):
    print(f"  Shape: {tensor.shape}")
    print(f"  NaN: {torch.isnan(tensor).any().item()}  Inf: {torch.isinf(tensor).any().item()}")
    print(f"  Max: {tensor.max().item()}  Min: {tensor.min().item()}  Mean: {tensor.mean().item()}  Std: {tensor.std().item()}")
