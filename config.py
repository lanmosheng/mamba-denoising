import yaml

import torch
from torchvision import transforms

import coder

# General config
# Datasets

def get_model(device, sampling_size=501):

    predictor = coder.MambaNet(num_patches=sampling_size)
    model = coder.Network(predictor, device=device)
    return model


