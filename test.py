import numpy as np

# 5001 方向，每个方向是一个 3D 向量，1000 个 patch
data = np.fromfile('dev/01t26t00.npy', dtype=np.float32).reshape((1000, 5001, 3))
label = np.fromfile('dev/F_01t26t00.npy', dtype=np.float32).reshape((1000, 3))

print(data.shape)   # (1000, 5001, 3)
print(label.shape)  # (1000, 3)
