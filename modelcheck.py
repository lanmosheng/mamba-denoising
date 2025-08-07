import torch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# ========== Load Data ==========
input_data = torch.load("debug_input.pt").cpu().float().numpy()   # [B, 1001, 3]
label_data = torch.load("debug_label.pt").cpu().float().numpy()   # [B, 3]

B, N, _ = input_data.shape
print(f"Loaded input_data shape: {input_data.shape}")
print(f"Loaded label_data shape: {label_data.shape}")

# ========== Analysis ==========
bad_samples = []

for i in range(B):
    patch = input_data[i]            # [1001, 3]
    label = label_data[i]            # [3]

    # 1. 均值 patch
    patch_mean = patch.mean(axis=0)  # [3]

    # 2. 单位化后计算夹角（cosine similarity）
    patch_mean_unit = patch_mean / (np.linalg.norm(patch_mean) + 1e-8)
    label_unit = label / (np.linalg.norm(label) + 1e-8)
    cos_sim = np.dot(patch_mean_unit, label_unit)
    angle = np.arccos(np.clip(cos_sim, -1.0, 1.0)) * 180 / np.pi

    # 3. patch 是否一致方向（方差极小）
    patch_std = np.std(patch, axis=0)
    patch_var = np.var(patch, axis=0)
    is_constant_patch = np.all(patch_var < 1e-6)

    # 4. patch 是否接近全零
    is_all_zero = np.all(np.abs(patch) < 1e-6)

    # 记录可疑样本
    if angle > 90 or is_constant_patch or is_all_zero:
        bad_samples.append((i, angle, is_constant_patch, is_all_zero))

    print(f"[{i:3}]  ∠(mean_patch, label) = {angle:.2f}° | Constant: {is_constant_patch} | Zero: {is_all_zero}")

print(f"\n⚠️ Found {len(bad_samples)} suspicious samples.\n")

# ========== Optional: Visualize a few bad samples ==========
VISUALIZE = True
SAVE_DIR = "bad_patch_vis"
os.makedirs(SAVE_DIR, exist_ok=True)

if VISUALIZE:
    for i, angle, constant, zero in bad_samples[:5]:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = input_data[i].T

        ax.quiver(np.zeros_like(x), np.zeros_like(y), np.zeros_like(z), x, y, z, length=0.1, normalize=True)
        ax.set_title(f"Sample {i} | ∠={angle:.1f}° | Const={constant} | Zero={zero}")
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.view_init(azim=60, elev=30)
        plt.savefig(f"{SAVE_DIR}/sample_{i}.png")
        plt.close()
