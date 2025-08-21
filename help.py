import os
import numpy as np

def get_total_cases(dataset_folder):
    total_cases = 0
    
    # 遍历 dataset 文件夹下的每个 mesh 文件夹
    for root, dirs, files in os.walk(dataset_folder):
        if 'lsd.npy' in files and 'gt.npy' in files:
            # 获取 lsd.npy 和 gt.npy 文件的路径
            lsd_path = os.path.join(root, 'lsd.npy')
            gt_path = os.path.join(root, 'gt.npy')
            
            # 读取 lsd 和 gt 数据
            train_data = np.load(lsd_path)  # shape: (nfaces, sampling_size, 3)
            train_label = np.load(gt_path)  # shape: (nfaces, 3)
            
            # 累加当前 mesh 文件夹中的样本数
            total_cases += train_data.shape[0]  # 每个 mesh 的样本数

    return total_cases

# 示例：调用函数并打印总的 case 数量
dataset_folder = 'train1601/ccylinder_n2'  # 替换为你实际的 dataset 路径
total_cases = get_total_cases(dataset_folder)
print(f'Total number of cases: {total_cases}')
