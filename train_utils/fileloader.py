# -*- coding: utf-8 -*-
import os
import json
from typing import Dict, List, Tuple, Optional
import numpy as np

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

class Loader:
    """
    Patch-centric loader.

    - 主数据是 patch 索引（patch_faces.npy），与 lsd/gt 分别位于独立目录：
        dataset_root/<mesh>/lsd.npy, gt.npy
        patch_root/<mesh>/patch_faces.npy  (K, M)  # 稀疏或稠密
    - sampling_size = lsd_r_size * lsd_t_size + 1
    - 每个 mesh 以“一个中心面=一个样本（其 patch）”展开
    - 支持 **稀疏 patch**：patch_faces.npy 只存被选中的 K 个中心 → 形状 (K, M)

    - 本 Loader **不做旋转与标准化**（已迁至 Trainer）。
    - 返回（默认）：
        data:  (num_batches, B, M, N, 3)
        label: (num_batches, B, M, 3)
      若 `return_face_idx=True`，额外返回：
        face_idx: (num_batches, B, M)  每个 (B,M) 位置对应的全局面号
    """
    def __init__(self,
                 dataset_root: str,
                 batch_size: int,
                 meta_path: Optional[str] = None,
                 patch_root: Optional[str] = None,
                 drop_last: bool = True,
                 shuffle_faces: bool = False,
                 mmap: bool = True,             # 兼容保留，Loader 内部不使用
                 return_face_idx: bool = False):   # ☆ 新增：是否返回 (B,M) 的面号
        self.dataset_root = dataset_root
        self.patch_root = patch_root
        if self.patch_root is None:
            raise ValueError("patch_root must be provided because patches and dataset live in separate folders")
        self.batch_size = int(batch_size)
        self.drop_last = bool(drop_last)
        self.shuffle_faces = bool(shuffle_faces)
        self.mmap_mode = 'r' if mmap else None

        # 新增
        self.return_face_idx = bool(return_face_idx)

        meta_path = meta_path or _default_meta_path(dataset_root)
        self.meta = _load_meta(meta_path)
        self.lsd_r_size = int(self.meta['lsd_r_size'])
        self.lsd_t_size = int(self.meta['lsd_t_size'])
        self.patch_num  = int(self.meta['patch_num'])
        self.sampling_size = self.lsd_r_size * self.lsd_t_size + 1

        self.mesh_dirs = _list_mesh_dirs(dataset_root)
        if not self.mesh_dirs:
            raise RuntimeError(f"No mesh dirs with lsd.npy & gt.npy under '{dataset_root}'")

        # Build records and validate shapes
        # record: (lsd_path, gt_path, patch_path, nfaces)
        self._records: List[Tuple[str, str, str, int]] = []
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
            # 允许 **稀疏** patch：只校验第二维等于 M
            if patches.ndim != 2 or patches.shape[1] != self.patch_num:
                raise ValueError(f"{patch_path}: expected (*,{self.patch_num}), got {patches.shape}")

            nfaces = int(lsd.shape[0])
            self._records.append((lsd_path, gt_path, patch_path, nfaces))

    def length(self) -> int:
        return len(self._records)

    def get_meta_data(self) -> Dict[str, int]:
        return {
            'lsd_r_size': self.lsd_r_size,
            'lsd_t_size': self.lsd_t_size,
            'patch_num': self.patch_num,
        }

    # ----------------------------
    # internal helpers
    # ----------------------------
    def _center_indices_sparse_or_dense(self, patch_faces: np.ndarray, nfaces: int, mesh_name: str) -> np.ndarray:
        """Return row indices to iterate patches.
        - 稠密模式：patch_faces.shape[0] == nfaces → 遍历所有行（等价于所有面）
        - 稀疏模式：patch_faces.shape[0] != nfaces → 每一行就是一个样本
        """
        ncenters = int(patch_faces.shape[0])
        dense = (ncenters == nfaces)
        if dense:
            center_idx = np.arange(nfaces, dtype=np.int64)
            print(f"[loader] mesh={mesh_name} patches=DENSE rows={ncenters} (== nfaces)", flush=True)
        else:
            center_idx = np.arange(ncenters, dtype=np.int64)
            print(f"[loader] mesh={mesh_name} patches=SPARSE rows={ncenters}/{nfaces}", flush=True)
        if self.shuffle_faces:
            np.random.shuffle(center_idx)
        return center_idx

    # ----------------------------
    # batch APIs
    # ----------------------------
    def generate_batch(self, mesh_idx: int, sampling_size: Optional[int] = None):
        """Load an entire mesh into numpy batches:
        Returns:
            data:  (num_batches, B, M, N, 3)
            label: (num_batches, B, M, 3)
            [face_idx: (num_batches, B, M)]  if return_face_idx=True
        """
        if sampling_size is not None and sampling_size != self.sampling_size:
            raise ValueError(f"sampling_size mismatch: got {sampling_size}, meta says {self.sampling_size}")

        lsd_path, gt_path, patch_path, nfaces = self._records[int(mesh_idx)]
        lsd         = np.load(lsd_path, mmap_mode=self.mmap_mode)
        gt          = np.load(gt_path,  mmap_mode=self.mmap_mode)
        patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode)

        if patch_faces.min() < 0 or patch_faces.max() >= nfaces:
            raise IndexError(f"Illegal face index in {patch_path}. Valid range [0,{nfaces-1}]")

        mesh_name = os.path.basename(os.path.dirname(patch_path))
        center_idx = self._center_indices_sparse_or_dense(patch_faces, nfaces, mesh_name)

        B = self.batch_size
        total = len(center_idx)
        num_batches = total // B
        if not self.drop_last and total % B != 0:
            num_batches += 1

        batches_X: List[np.ndarray] = []
        batches_Y: List[np.ndarray] = []
        batches_I: List[np.ndarray] = [] if self.return_face_idx else None

        for bi in range(num_batches):
            start = bi * B
            end = min(start + B, total)
            cur_centers = center_idx[start:end]
            cur_B = len(cur_centers)
            if cur_B < B and self.drop_last:
                break

            Xb = np.empty((cur_B, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            Yb = np.empty((cur_B, self.patch_num, 3), dtype=np.float32)
            if self.return_face_idx:
                Ib = np.empty((cur_B, self.patch_num), dtype=np.int64)

            for i, c in enumerate(cur_centers):
                faces = patch_faces[c]            # (M,)
                X = lsd[faces]                    # (M,N,3)
                Y = gt[faces]                     # (M,3)

                Xb[i] = X.astype(np.float32)
                Yb[i] = Y.astype(np.float32)
                if self.return_face_idx:
                    Ib[i] = faces.astype(np.int64)

            batches_X.append(Xb)
            batches_Y.append(Yb)
            if self.return_face_idx:
                batches_I.append(Ib)

        if not batches_X:
            empty_X = np.empty((0, self.batch_size, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            empty_Y = np.empty((0, self.batch_size, self.patch_num, 3), dtype=np.float32)
            if self.return_face_idx:
                empty_I = np.empty((0, self.batch_size, self.patch_num), dtype=np.int64)
                return empty_X, empty_Y, empty_I
            return empty_X, empty_Y

        data  = np.stack(batches_X, axis=0)
        label = np.stack(batches_Y, axis=0)
        if self.return_face_idx:
            face_idx = np.stack(batches_I, axis=0)
            return data, label, face_idx
        return data, label

    def iter_batches(self, mesh_idx: int, sampling_size: int | None = None):
        """Yield small batches for a mesh (streaming).
        Yields:
            Xb: (B, M, N, 3)
            Yb: (B, M, 3)
            [Ib: (B, M)] if return_face_idx=True
        """
        if sampling_size is not None and sampling_size != self.sampling_size:
            raise ValueError(f"sampling_size mismatch: got {sampling_size}, meta says {self.sampling_size}")

        lsd_path, gt_path, patch_path, nfaces = self._records[int(mesh_idx)]
        lsd         = np.load(lsd_path, mmap_mode=self.mmap_mode)
        gt          = np.load(gt_path,  mmap_mode=self.mmap_mode)
        patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode)

        if patch_faces.min() < 0 or patch_faces.max() >= nfaces:
            raise IndexError(f"Illegal face index in {patch_path}. Valid range [0,{nfaces-1}]")

        mesh_name = os.path.basename(os.path.dirname(patch_path))
        center_idx = self._center_indices_sparse_or_dense(patch_faces, nfaces, mesh_name)

        B = self.batch_size
        total = len(center_idx)
        num_batches = total // B + (0 if self.drop_last or total % B == 0 else 1)

        for bi in range(num_batches):
            start = bi * B
            end = min(start + B, total)
            if end - start < B and self.drop_last:
                break
            cur_centers = center_idx[start:end]
            cur_B = len(cur_centers)

            Xb = np.empty((cur_B, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            Yb = np.empty((cur_B, self.patch_num, 3), dtype=np.float32)
            if self.return_face_idx:
                Ib = np.empty((cur_B, self.patch_num), dtype=np.int64)

            for i, c in enumerate(cur_centers):
                faces = patch_faces[c]                # (M,)
                X = lsd[faces]                        # (M,N,3)
                Y = gt[faces]                         # (M,3)

                Xb[i] = X.astype(np.float32)
                Yb[i] = Y.astype(np.float32)
                if self.return_face_idx:
                    Ib[i] = faces.astype(np.int64)

            if self.return_face_idx:
                yield Xb, Yb, Ib
            else:
                yield Xb, Yb

    def count_batches(self, mesh_idx: int) -> int:
        """Return number of batches for this mesh, honoring sparse patch_faces and drop_last."""
        _, _, patch_path, nfaces = self._records[int(mesh_idx)]
        patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode)
        ncenters = int(patch_faces.shape[0])
        total = ncenters  # 稀疏/稠密统一按行数
        B = self.batch_size
        return total // B + (0 if self.drop_last or total % B == 0 else 1)
