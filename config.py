import yaml

import torch
from torchvision import transforms

import coder

# General config
# Datasets

def get_model(device, sampling_size=1601, lsd_r = 40, lsd_t = 40):

    predictor = coder.MambaNet(num_patches=sampling_size,lsd_r_size=lsd_r, lsd_t_size=lsd_t, add_pos_emb=False)
    model = coder.Network(predictor, device=device)
    return model


