import os
import numpy as np

# 替换为你的目录路径
data_dir = "/home/jzw/Proj/mamba_denoising/train"

# 遍历所有 .npy 文件
for fname in os.listdir(data_dir):
    if not fname.endswith(".npy"):
        continue
    if fname.startswith("F"):  # 跳过 F 文件
        continue
    fpath = os.path.join(data_dir, fname)
    try:
        arr = np.fromfile(fpath, dtype=np.float32)
        if arr.size % (1001 * 3) != 0:
            print(f"[❌] {fname}: {arr.size} elements (NOT divisible by 1001*3)")
        else:
            reshaped = arr.reshape((-1, 1001, 3))
            print(f"[✅] {fname}: OK, shape = {reshaped.shape}")
    except Exception as e:
        print(f"[‼️] Error reading {fname}: {e}")
