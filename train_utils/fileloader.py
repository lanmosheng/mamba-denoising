# -*- coding: utf-8 -*-
import os
import json
from typing import Dict, List, Tuple, Optional
import numpy as np

# ----------------------------
# Utility functions
# ----------------------------

def _list_mesh_dirs(dataset_root: str) -> List[str]:
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


def _safe_norm(v: np.ndarray, eps: float = 1e-8) -> float:
    """Return max(||v||, eps) as a *lower-bounded* norm for stable normalization."""
    n = float(np.linalg.norm(v))
    return n if n > eps else eps


def _rotation_matrix_from_a_to_b(a: np.ndarray, b: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Return 3x3 rotation matrix that rotates vector a to b."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    an = a / _safe_norm(a, eps)
    bn = b / _safe_norm(b, eps)

    v = np.cross(an, bn)       # rotation axis * sin(theta)
    c = float(np.dot(an, bn))  # cos(theta)
    s = float(np.linalg.norm(v))  # |sin(theta)| -- use *true* norm for branch

    if s < 1e-12:
        if c > 0.0:  # parallel
            return np.eye(3, dtype=np.float64)
        # opposite: 180° rotate around any orthogonal axis
        axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(an[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        axis = axis / _safe_norm(axis, eps)
        # Rodrigues for 180°: R = I + 2*K^2 with s=0, c=-1  (since sin=0, 1-c = 2)
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]], dtype=np.float64)
        return np.eye(3, dtype=np.float64) + 2.0 * (K @ K)

    axis = v / s
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]], dtype=np.float64)
    # Rodrigues: R = I + K*sin + K^2*(1-cos)
    R = np.eye(3, dtype=np.float64) + K * s + (K @ K) * (1.0 - c)
    return R


def _zscore_patch(X: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Patch-level z-score over (M*N) jointly per channel of last dim=3. X: (M, N, 3)"""
    flat = X.reshape(-1, 3)
    mu = flat.mean(axis=0)
    std = flat.std(axis=0)
    std = np.maximum(std, eps)
    return (X - mu) / std

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
    - generate_batch 返回：
        data:  (num_batches, B, M, N, 3)
        label: (num_batches, B, M, 3)
    - 返回前做 patch 级 z-score + patch-wise 旋转到 (1,0,0)
    """
    def __init__(self,
                 dataset_root: str,
                 batch_size: int,
                 meta_path: Optional[str] = None,
                 patch_root: Optional[str] = None,
                 drop_last: bool = True,
                 shuffle_faces: bool = False,
                 mmap: bool = True,
                 rotation_anchor: str = 'center',  # 'center' or 'mean'
                 eps: float = 1e-6):
        self.dataset_root = dataset_root
        self.patch_root = patch_root
        if self.patch_root is None:
            raise ValueError("patch_root must be provided because patches and dataset live in separate folders")
        self.batch_size = int(batch_size)
        self.drop_last = bool(drop_last)
        self.shuffle_faces = bool(shuffle_faces)
        self.mmap_mode = 'r' if mmap else None
        self.rotation_anchor = rotation_anchor
        self.eps = eps

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
        """Load an entire mesh into numpy batches: returns (data_batches, label_batches)."""
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

        target = np.array([1.0, 0.0, 0.0], dtype=np.float64)  # +X

        for bi in range(num_batches):
            start = bi * B
            end = min(start + B, total)
            cur_centers = center_idx[start:end]
            cur_B = len(cur_centers)
            if cur_B < B and self.drop_last:
                break

            Xb = np.empty((cur_B, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            Yb = np.empty((cur_B, self.patch_num, 3), dtype=np.float32)

            for i, c in enumerate(cur_centers):
                faces = patch_faces[c]            # (M,)
                X = lsd[faces]                    # (M,N,3)
                Y = gt[faces]                     # (M,3)

                # anchor selection from observed LSD (no GT)
                center_dirs = X[:, self.sampling_size - 1]
                if self.rotation_anchor == 'center':
                    anchor = center_dirs[0]
                else:  # 'mean'
                    center_norm = center_dirs / np.maximum(
                        np.linalg.norm(center_dirs, axis=1, keepdims=True), self.eps
                    )
                    anchor = center_norm.mean(axis=0)

                # degenerate fallback: zero/NaN anchor -> use mean valid; still bad -> +X
                if np.linalg.norm(anchor) < 1e-6 or not np.isfinite(anchor).all():
                    mask = np.isfinite(center_dirs).all(axis=1) & (
                        np.linalg.norm(center_dirs, axis=1) > 1e-6
                    )
                    anchor = (
                        center_dirs[mask].mean(axis=0)
                        if mask.any()
                        else np.array([1.0, 0.0, 0.0], dtype=np.float64)
                    )

                R = _rotation_matrix_from_a_to_b(anchor, target, eps=self.eps)

                # rotate
                X_rot = X.reshape(-1, 3) @ R.T
                X_rot = X_rot.reshape(self.patch_num, self.sampling_size, 3)
                Y_rot = Y @ R.T

                # patch-level z-score
                Xn = _zscore_patch(X_rot, eps=self.eps)

                Xb[i] = Xn.astype(np.float32)
                Yb[i] = Y_rot.astype(np.float32)

            batches_X.append(Xb)
            batches_Y.append(Yb)

        if not batches_X:
            return (
                np.empty((0, self.batch_size, self.patch_num, self.sampling_size, 3), dtype=np.float32),
                np.empty((0, self.batch_size, self.patch_num, 3), dtype=np.float32),
            )

        data  = np.stack(batches_X, axis=0)
        label = np.stack(batches_Y, axis=0)
        return data, label

    def iter_batches(self, mesh_idx: int, sampling_size: int | None = None):
        """Yield small batches for a mesh (streaming). Shape same as generate_batch per item."""
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

        target = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        for bi in range(num_batches):
            start = bi * B
            end = min(start + B, total)
            if end - start < B and self.drop_last:
                break
            cur_centers = center_idx[start:end]
            cur_B = len(cur_centers)

            Xb = np.empty((cur_B, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            Yb = np.empty((cur_B, self.patch_num, 3), dtype=np.float32)

            for i, c in enumerate(cur_centers):
                faces = patch_faces[c]                # (M,)
                X = lsd[faces]                        # (M,N,3)
                Y = gt[faces]                         # (M,3)

                # anchor selection from observed LSD (no GT)
                center_dirs = X[:, self.sampling_size - 1]
                if self.rotation_anchor == 'center':
                    anchor = center_dirs[0]
                else:
                    center_norm = center_dirs / np.maximum(
                        np.linalg.norm(center_dirs, axis=1, keepdims=True), self.eps
                    )
                    anchor = center_norm.mean(axis=0)

                if np.linalg.norm(anchor) < 1e-6 or not np.isfinite(anchor).all():
                    mask = np.isfinite(center_dirs).all(axis=1) & (
                        np.linalg.norm(center_dirs, axis=1) > 1e-6
                    )
                    anchor = (
                        center_dirs[mask].mean(axis=0)
                        if mask.any()
                        else np.array([1.0, 0.0, 0.0], dtype=np.float64)
                    )

                R = _rotation_matrix_from_a_to_b(anchor, target, eps=self.eps)

                X_rot = X.reshape(-1, 3) @ R.T
                X_rot = X_rot.reshape(self.patch_num, self.sampling_size, 3)
                Y_rot = Y @ R.T

                Xn = _zscore_patch(X_rot, eps=self.eps)

                Xb[i] = Xn.astype(np.float32)
                Yb[i] = Y_rot.astype(np.float32)

            # stream out
            yield Xb, Yb

    def count_batches(self, mesh_idx: int) -> int:
        """Return number of batches for this mesh, honoring sparse patch_faces and drop_last."""
        _, _, patch_path, nfaces = self._records[int(mesh_idx)]
        patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode)
        ncenters = int(patch_faces.shape[0])
        total = ncenters  # 稀疏/稠密统一按行数
        B = self.batch_size
        return total // B + (0 if self.drop_last or total % B == 0 else 1)
