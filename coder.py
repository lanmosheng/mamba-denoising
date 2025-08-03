import torch.nn as nn
import torch
import torch.nn.functional as F
import numpy as np
from mamba_ssm import Mamba

class Network(nn.Module):

    def __init__(self, predictor, device):
        super().__init__()

        self.predictor = predictor.to(device)
        self._device = device



    def compute_loss(self, data, lab):
        output = F.normalize(self.pred(data), dim=1)
        loss = torch.nn.functional.cosine_similarity(output, lab, dim=1)
        loss2=torch.acos(loss)*180.0/np.pi
        isnan=torch.isnan(loss2)
        for i in range(isnan.shape[0]):
            if isnan[i]==True:
                loss2[i]=0
        return loss2

    def pred(self, data):
        output = self.predictor(data)
        return output


class MambaNet(nn.Module):
    def __init__(self, num_patches=1001, d_model=64, d_state=16, d_conv=4, expand=2, depth=6):
        super().__init__()
        self.num_patches = num_patches
        self.input_proj = nn.Linear(3, d_model)

        self.mamba_layers = nn.ModuleList([
            Mamba(d_model=d_model, d_state=d_state, d_conv=d_conv, expand=expand)
            for _ in range(depth)
        ])

        self.output_head = nn.Sequential(
            nn.Linear(d_model, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Linear(128, 3),
            nn.Tanh()
        )

    def forward(self, x):
        # x: [B, num_patches, 3]
        assert x.shape[1] == self.num_patches, \
            f"Expected {self.num_patches} patches, but got {x.shape[1]}"

        x = self.input_proj(x)  # [B, num_patches, d_model]

        for layer in self.mamba_layers:
            x = layer(x)

        x = x.mean(dim=1)  # 聚合每个 patch 表征 → [B, d_model]
        x = self.output_head(x)  # → [B, 3]
        return x