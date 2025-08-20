import os
import numpy as np
import json

def read_filename(data_path):
    mesh_folders = []
    for root, dirs, files in os.walk(data_path):
        # 寻找每个 mesh 文件夹中包含的 lsd.npy 和 gt.npy
        if 'lsd.npy' in files and 'gt.npy' in files:
            mesh_folders.append(root)
    return mesh_folders

def load_meta(meta_path):
    """加载 meta.json 文件并返回字典"""
    with open(meta_path, 'r') as f:
        meta = json.load(f)
    return meta

class Loader():
    def __init__(self, dataset_folder, batchsize, meta_path):
        # 读取 mesh 文件夹
        mesh_folders = read_filename(dataset_folder)

        # 加载 meta 数据
        self.meta = load_meta(meta_path)
        
        # 从 meta 数据中提取 lsd_r_size, lsd_t_size, patch_num
        self.lsd_r_size = self.meta["lsd_r_size"]
        self.lsd_t_size = self.meta["lsd_t_size"]
        self.patch_num = self.meta["patch_num"]

        # 计算 sampling_size
        self.sampling_size = self.lsd_r_size * self.lsd_t_size + 1
        print(f"Calculated sampling_size: {self.sampling_size}")

        # 存储每个 mesh 文件夹的 lsd.npy 和 gt.npy
        self.mesh_data = []
        for folder in mesh_folders:
            lsd_path = os.path.join(folder, 'lsd.npy')
            gt_path = os.path.join(folder, 'gt.npy')
            self.mesh_data.append((lsd_path, gt_path))
        
        self.batchsize = batchsize
        self.train_num = len(self.mesh_data)
        
    def length(self):
        return self.train_num
    
    def generate_batch(self, idx):
        # 获取当前 mesh 文件夹中的 lsd 和 gt 数据
        lsd_path, gt_path = self.mesh_data[idx]
        
        # 读取 lsd 和 gt 数据
        train_data = np.load(lsd_path)  # shape: (nfaces, sampling_size, 3)
        train_label = np.load(gt_path)  # shape: (nfaces, 3)
        
        # 校验数据的形状是否正确
        if train_data.shape[1] != self.sampling_size or train_data.shape[2] != 3:
            raise ValueError(f"Shape mismatch in {lsd_path}: expected shape (nfaces, {self.sampling_size}, 3), got {train_data.shape}")
        
        # 获取数据样本数量
        case_number = train_data.shape[0]
        
        # 计算有效的批次数量，舍弃最后一个不满的批次
        batch_number = case_number // self.batchsize
        
        # 截取有效的训练数据和标签
        train_data = train_data[:batch_number * self.batchsize].reshape((batch_number, self.batchsize, self.sampling_size, 3))
        train_label = train_label[:batch_number * self.batchsize].reshape((batch_number, self.batchsize, 3))
        
        return train_data, train_label

    def get_meta_data(self):
        """返回 meta 数据"""
        return {
            'lsd_r_size': self.lsd_r_size,
            'lsd_t_size': self.lsd_t_size,
            'patch_num': self.patch_num
        }
