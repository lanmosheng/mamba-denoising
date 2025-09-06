# -*- coding: utf-8 -*-
import os
import json
from typing import Dict, List, Tuple, Optional
import numpy as np
import math
import random
# ----------------------------
# Utility functions
# ----------------------------

def _list_mesh_dirs(dataset_root: str) -> List[str]:
    """Recursively collect subdirs that contain both lsd.npy and gt.npy."""
    mesh_dirs = []
    for root, dirs, files in os.walk(dataset_root):
        if 'lsd.npy' in files and 'gt.npy' in files:
            mesh_dirs.append(root)
    mesh_dirs.sort()
    return mesh_dirs


def _default_meta_path(dataset_root: str) -> str:
    return os.path.join(dataset_root, 'meta.json')


def _load_meta(meta_path: str) -> Dict:
    with open(meta_path, 'r') as f:
        meta = json.load(f)
    for k in ['lsd_r_size', 'lsd_t_size', 'patch_num']:
        if k not in meta:
            raise KeyError(f"Missing '{k}' in meta.json")
    return meta


def _resolve_patch_path(mesh_dir: str, patch_root: Optional[str] = None) -> str:
    """
    Resolve patch_faces.npy path for a given mesh:
      REQUIRED: patch_root must be provided and point to a separate folder tree,
      with layout: <patch_root>/<mesh_name>/patch_faces.npy
    """
    if patch_root is None:
        raise ValueError("patch_root is required (patches and dataset are in separate folders)")
    mesh_name = os.path.basename(mesh_dir.rstrip(os.sep))
    cand = os.path.join(patch_root, mesh_name, 'patch_faces.npy')
    if os.path.exists(cand):
        return cand
    raise FileNotFoundError(f"patch_faces.npy not found for mesh '{mesh_name}' under patch_root='{patch_root}'")


# ----------------------------
# Loader
# ----------------------------
import os, math
from typing import Optional, List, Tuple
import numpy as np

class S1FaceLoaderV2:
    """
    第一层 FaceEncoder 训练数据管道（仅旋转规范化；不做 z-score）
    - 旧接口保留：dataset_root, batch_size, meta_path=None, patch_root=None
    - 文件搜索/校验沿用：_default_meta_path/_load_meta/_list_mesh_dirs/_resolve_patch_path
    - 训练逻辑：
        1) 每个 mesh 内，打乱所有 patch 行顺序
        2) 每个 epoch 仅取该 mesh 的一个切片（如 10%），跨 epoch 不重叠；确保每个 epoch 每个 mesh 至少 1 个 patch
        3) 以 patch 为单位：取中心面（faces[0]）的 LSD token0 作为朝向锚点，计算 R 使其对齐到 (1,0,0)
        4) 用同一 R 旋转该 patch 的所有 LSD；并用同一 R 旋转对应的 GT 法向
        5) 逐面打包为 batch: X ∈ [B, N, 3], Y ∈ [B, 3]
    """
    def __init__(
        self,
        dataset_root: str,
        batch_size: int,
        meta_path: Optional[str] = None,
        patch_root: Optional[str] = None,
        *,
        slice_ratio: float = 0.10,         # 每个 epoch 取多少比例的 patch（默认 10%）
        shuffle_faces: bool = True,
        shuffle_seed: Optional[int] = 0,
        drop_last: bool = True,
        mmap_mode: Optional[str] = 'r',
    ):
        # ---- 旧工具函数：meta & 路径解析 ----
        meta_path = meta_path or _default_meta_path(dataset_root)
        self.meta = _load_meta(meta_path)
        self.lsd_r_size   = int(self.meta['lsd_r_size'])
        self.lsd_t_size   = int(self.meta['lsd_t_size'])
        self.patch_num    = int(self.meta['patch_num'])
        self.sampling_size = 1 + self.lsd_r_size * self.lsd_t_size  # N

        self.dataset_root = dataset_root
        self.patch_root   = patch_root
        self.batch_size   = int(batch_size)
        self.slice_ratio  = float(slice_ratio)
        assert 0 < self.slice_ratio <= 1.0, "slice_ratio 必须在 (0,1] 内"
        self.shuffle_faces = bool(shuffle_faces)
        self.shuffle_seed  = shuffle_seed
        self.drop_last     = bool(drop_last)
        self.mmap_mode     = mmap_mode

        # ---- 旧文件搜索逻辑 ----
        self.mesh_dirs = _list_mesh_dirs(dataset_root)
        if not self.mesh_dirs:
            raise RuntimeError(f"No mesh dirs with lsd.npy & gt.npy under '{dataset_root}'")

        self._records: List[Tuple[str, str, str, int]] = []  # (lsd_path, gt_path, patch_path, nfaces)

        for mdir in self.mesh_dirs:
            lsd_path = os.path.join(mdir, 'lsd.npy')
            gt_path  = os.path.join(mdir, 'gt.npy')
            patch_path = _resolve_patch_path(mdir, patch_root=self.patch_root)

            lsd     = np.load(lsd_path, mmap_mode=self.mmap_mode)
            gt      = np.load(gt_path,  mmap_mode=self.mmap_mode)
            patches = np.load(patch_path, mmap_mode=self.mmap_mode)

            if lsd.ndim != 3 or lsd.shape[1] != self.sampling_size or lsd.shape[2] != 3:
                raise ValueError(f"{lsd_path}: expected (nfaces,{self.sampling_size},3), got {lsd.shape}")
            if gt.ndim != 2 or gt.shape[1] != 3:
                raise ValueError(f"{gt_path}: expected (nfaces,3), got {gt.shape}")
            if patches.ndim != 2 or patches.shape[1] != self.patch_num:
                raise ValueError(f"{patch_path}: expected (K,{self.patch_num}), got {patches.shape}")

            nfaces = int(lsd.shape[0])
            self._records.append((lsd_path, gt_path, patch_path, nfaces))

        self._mesh_order = list(range(len(self._records)))
        # ---- 每个 mesh 的打乱顺序与切片大小 ----
        self._mesh_patch_orders: List[np.ndarray] = []  # 打乱后的 patch 行号
        self._mesh_k: List[int] = []                   # 每个 mesh 的 patch 行数 K
        self._mesh_slice_size: List[int] = []          # 每片大小（>=1）
        self._mesh_slices_per_mesh: List[int] = []     # 每个 mesh 的切片数（>=1）

        for mi, (_, _, patch_path, _) in enumerate(self._records):
            pf = np.load(patch_path, mmap_mode=self.mmap_mode)
            k = int(pf.shape[0])
            order = np.arange(k, dtype=np.int64)
            if self.shuffle_faces:
                rng = np.random.RandomState(None if self.shuffle_seed is None else (self.shuffle_seed + mi))
                rng.shuffle(order)
            self._mesh_patch_orders.append(order)
            self._mesh_k.append(k)
            slice_size = max(1, int(math.ceil(k * self.slice_ratio)))
            self._mesh_slice_size.append(slice_size)
            self._mesh_slices_per_mesh.append(max(1, int(math.ceil(k / slice_size))))

        self._epoch = 0

        # 旋转目标方向（固定到 (1,0,0)）
        self._canonical_dir = np.array([1.0, 0.0, 0.0], dtype=np.float64)

    # ========== 公共接口 ==========
    def set_epoch(self, epoch: int):
        self._epoch = int(epoch)
        random.shuffle(self._mesh_order)

    def iter_batches(self):
        """
        逐批产出单面样本（已旋转规范化）：
          X: (B, N, 3)   —— 旋转后的 LSD
          Y: (B, 3)      —— 旋转后的 GT 法向
          元信息：mesh_idx / face_idx / center_face / patch_row
        """
        buf_X, buf_Y = [], []
        buf_mid, buf_fid, buf_cen, buf_row = [], [], [], []

        for mi in self._mesh_order:
            lsd_path, gt_path, patch_path, nfaces = self._records[mi]
            lsd         = np.load(lsd_path, mmap_mode=self.mmap_mode)  # (nfaces, N, 3)
            gt          = np.load(gt_path,  mmap_mode=self.mmap_mode)  # (nfaces, 3)
            patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode) # (K, M)

            if patch_faces.min() < 0 or patch_faces.max() >= nfaces:
                raise IndexError(f"Illegal face index in {patch_path}. Valid range [0,{nfaces-1}]")

            center_rows = self._center_indices_sparse_or_dense(
                patch_faces, nfaces, os.path.basename(os.path.dirname(patch_path))
            )

            k = self._mesh_k[mi]
            slice_size       = self._mesh_slice_size[mi]
            slices_per_mesh  = self._mesh_slices_per_mesh[mi]
            sid = self._epoch % slices_per_mesh
            beg = sid * slice_size
            end = min(k, beg + slice_size)
            chosen_rows = center_rows[beg:end]  # 保证至少 1 行

            for r in chosen_rows:
                faces = patch_faces[int(r)]          # (M,)
                # 以中心面（faces[0]）LSD 的 token0 作为锚点，估计旋转
                X_patch = lsd[faces]                 # (M, N, 3)
                R = self._rotation_from_center_token0(X_patch)  # (3,3)
                X_rot = X_patch @ R.T                # 旋转整个 patch 的 LSD -> (M, N, 3)

                # 标签也用同一 R 旋转到同一坐标系
                Y_rot = gt[faces] @ R.T              # (M, 3)

                center = int(faces[0])
                for i, fid in enumerate(faces):
                    buf_X.append(X_rot[i].astype(np.float32))
                    buf_Y.append(Y_rot[i].astype(np.float32))
                    buf_mid.append(mi)
                    buf_fid.append(int(fid))
                    buf_cen.append(center)
                    buf_row.append(int(r))

                    if len(buf_X) == self.batch_size:
                        yield self._pack_batch(buf_X, buf_Y, buf_mid, buf_fid, buf_cen, buf_row)
                        buf_X, buf_Y, buf_mid, buf_fid, buf_cen, buf_row = [], [], [], [], [], []

        if len(buf_X) and not self.drop_last:
            yield self._pack_batch(buf_X, buf_Y, buf_mid, buf_fid, buf_cen, buf_row)

    # ========== 内部工具 ==========
    def _pack_batch(self, X_list, Y_list, mid_list, fid_list, cen_list, row_list):
        X = np.stack(X_list, axis=0)  # (B, N, 3) float32
        Y = np.stack(Y_list, axis=0)  # (B, 3)    float32
        return {
            'X': X, 'Y': Y,
            'mesh_idx': np.asarray(mid_list, dtype=np.int64),
            'face_idx': np.asarray(fid_list, dtype=np.int64),
            'center_face': np.asarray(cen_list, dtype=np.int64),
            'patch_row': np.asarray(row_list, dtype=np.int64),
        }

    def _center_indices_sparse_or_dense(self, patch_faces: np.ndarray, nfaces: int, mesh_name: str) -> np.ndarray:
        ncenters = int(patch_faces.shape[0])
        dense = (ncenters == nfaces)
        if dense:
            idx = np.arange(nfaces, dtype=np.int64)
            print(f"[loader] mesh={mesh_name} patches=DENSE rows={nfaces}", flush=True)
        else:
            idx = np.arange(ncenters, dtype=np.int64)
            print(f"[loader] mesh={mesh_name} patches=SPARSE rows={ncenters}/{nfaces}", flush=True)
        if self.shuffle_faces:
            np.random.shuffle(idx)
        return idx

    # === 旋转估计 ===
    def _rotation_from_center_token0(self, X_patch: np.ndarray) -> np.ndarray:
        """
        用中心面（faces[0]）的 LSD token0 作为朝向向量 a，计算 R 使 a -> (1,0,0)。
        仅依赖 LSD，不触及 GT。对退化情形做稳定处理。
        """
        a = X_patch[0, 0].astype(np.float64)           # 中心面的 token0
        a_norm = np.linalg.norm(a)
        if a_norm < 1e-12:
            return np.eye(3, dtype=np.float64)         # 退化：不旋转
        a /= a_norm
        b = self._canonical_dir                        # (1,0,0)
        return self._rotation_matrix_a2b(a, b)

    @staticmethod
    def _rotation_matrix_a2b(a: np.ndarray, b: np.ndarray, eps: float = 1e-8) -> np.ndarray:
        """返回将单位向量 a 旋到 b 的 3x3 旋转矩阵（稳定处理平行/反平行）。"""
        v = np.cross(a, b)
        c = float(np.dot(a, b))    # cosθ
        s = np.linalg.norm(v)      # |v| = sinθ

        if s < eps:
            if c > 0.0:
                return np.eye(3, dtype=np.float64)  # a≈b
            # a≈-b：选任一与 a 不共线的轴做 180°
            axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
            if abs(a[0]) > 0.9:
                axis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
            axis = axis - a * np.dot(a, axis)
            axis /= (np.linalg.norm(axis) + eps)
            K = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]], dtype=np.float64)
            return np.eye(3) + 2 * (K @ K)

        vx, vy, vz = v / s
        K = np.array([[0, -vz, vy],
                      [vz, 0, -vx],
                      [-vy, vx, 0]], dtype=np.float64)
        return np.eye(3) + K * s + (K @ K) * ((1.0 - c) / (s * s + eps))
    def total_batches(self) -> int:
        """
        返回“当前 epoch（已 set_epoch）”会产出的 batch 总数。
        计算基于：每个 mesh 的切片大小 × patch_num（每行patch产出M个面样本） × batch_size。
        """
        total_samples = 0
        for mi in range(len(self._records)):
            k = self._mesh_k[mi]
            slice_size = self._mesh_slice_size[mi]
            slices_per_mesh = self._mesh_slices_per_mesh[mi]
            sid = self._epoch % slices_per_mesh
            beg = sid * slice_size
            end = min(k, beg + slice_size)
            chosen_rows = max(0, end - beg)
            total_samples += chosen_rows * self.patch_num  # 每行patch产出 M 个“面”样本

        bs = self.batch_size
        if self.drop_last:
            return int(total_samples // bs)
        else:
            return int((total_samples + bs - 1) // bs)