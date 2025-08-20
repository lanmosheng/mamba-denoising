import os
import numpy as np

def read_filename(data_path):
    mesh_folders = []
    for root, dirs, files in os.walk(data_path):
        # 寻找每个 mesh 文件夹中包含的 lsd.npy 和 gt.npy
        if 'lsd.npy' in files and 'gt.npy' in files:
            mesh_folders.append(root)
    return mesh_folders

class Loader():
    def __init__(self, dataset_folder, batchsize):
        # 读取 mesh 文件夹
        mesh_folders = read_filename(dataset_folder)

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
    
    def generate_batch(self, idx, sampling_size: int):
        # 获取当前 mesh 文件夹中的 lsd 和 gt 数据
        lsd_path, gt_path = self.mesh_data[idx]
        
        # 读取 lsd 和 gt 数据
        train_data = np.load(lsd_path)  # shape: (nfaces, sampling_size, 3)
        train_label = np.load(gt_path)  # shape: (nfaces, 3)
        
        # 获取数据样本数量
        case_number = train_data.shape[0]
        
        # 计算有效的批次数量，舍弃最后一个不满的批次
        batch_number = case_number // self.batchsize
        
        # 截取有效的训练数据和标签
        train_data = train_data[:batch_number * self.batchsize].reshape((batch_number, self.batchsize, sampling_size, 3))
        train_label = train_label[:batch_number * self.batchsize].reshape((batch_number, self.batchsize, 3))
        
        return train_data, train_label

