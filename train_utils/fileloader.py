import os
import json
import math
from typing import Dict, List, Tuple, Optional

import numpy as np


def _list_mesh_dirs(dataset_root: str) -> List[str]:
    """Return mesh folders that contain lsd.npy and gt.npy."""
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
    # basic checks
    for k in ['lsd_r_size', 'lsd_t_size', 'patch_num']:
        if k not in meta:
            raise KeyError(f"Missing '{k}' in meta.json")
    return meta


def _resolve_patch_path(mesh_dir: str, patch_root: Optional[str] = None, dataset_root: Optional[str] = None) -> str:
    """Resolve path to patch_faces.npy for a given mesh directory.
    Priority:
      1) <patch_root>/<mesh_name>/patch_faces.npy  (if patch_root provided)
      2) <dataset_root>/patches/<mesh_name>/patch_faces.npy (if dataset_root provided)
      3) <mesh_dir>/patches/patch_faces.npy
    """
    mesh_name = os.path.basename(mesh_dir.rstrip(os.sep))
    tried = []

    if patch_root is not None:
        cand = os.path.join(patch_root, mesh_name, 'patch_faces.npy')
        tried.append(cand)
        if os.path.exists(cand):
            return cand

    if dataset_root is not None:
        cand = os.path.join(dataset_root, 'patches', mesh_name, 'patch_faces.npy')
        tried.append(cand)
        if os.path.exists(cand):
            return cand

    cand = os.path.join(mesh_dir, 'patches', 'patch_faces.npy')
    tried.append(cand)
    if os.path.exists(cand):
        return cand

    tried_list = "\n- ".join(tried)
    raise FileNotFoundError(
        f"patch_faces.npy not found for mesh '{mesh_name}'. Tried:\n- {tried_list}")


def _safe_norm(v: np.ndarray, eps: float = 1e-8) -> float:
    return float(np.linalg.norm(v) + eps)


def _rotation_matrix_from_a_to_b(a: np.ndarray, b: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Return 3x3 rotation matrix that rotates vector a to vector b.
    Handles near-parallel and near-opposite cases robustly.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    an = a / _safe_norm(a, eps)
    bn = b / _safe_norm(b, eps)

    v = np.cross(an, bn)
    c = float(np.dot(an, bn))  # cos(theta)
    s = _safe_norm(v, eps)      # |sin(theta)|

    if s < 1e-12:
        # vectors are parallel or anti-parallel
        if c > 0.0:  # nearly identical
            return np.eye(3, dtype=np.float64)
        # opposite direction: rotate 180° around any axis orthogonal to 'an'
        # find an orthogonal axis
        axis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        if abs(an[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        u = axis - an * float(np.dot(axis, an))
        u = u / _safe_norm(u, eps)
        # Rodrigues for 180°: R = I + 2[K(u)]^2 where K is cross-product matrix
        ux, uy, uz = u
        K = np.array([[0, -uz, uy], [uz, 0, -ux], [-uy, ux, 0]], dtype=np.float64)
        return (np.eye(3, dtype=np.float64) + 2.0 * (K @ K))

    # Rodrigues' rotation formula
    vx, vy, vz = v / s
    K = np.array([[0, -vz, vy], [vz, 0, -vx], [-vy, vx, 0]], dtype=np.float64)
    R = np.eye(3, dtype=np.float64) + K * s + (K @ K) * (1.0 - c)
    return R


def _zscore_patch(X: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Patch-level z-score over (M*N) jointly per channel of last dim=3.
    X: (M, N, 3)
    """
    flat = X.reshape(-1, 3)
    mu = flat.mean(axis=0)
    std = flat.std(axis=0)
    std = np.maximum(std, eps)
    return (X - mu) / std


class Loader:
    """
    Patch-centric loader.

    - Primary data is patch indices (patch_faces.npy). LSD/GT are lazily sliced.
    - Computes sampling_size from meta: lsd_r_size * lsd_t_size + 1.
    - For each mesh, treats every face as one training sample (its patch).
    - Returns per-mesh batches shaped:
        data:  (num_batches, batch_size, M, N, 3)
        label: (num_batches, batch_size, M, 3)
    - Applies patch-level z-score and patch-wise rotation-to-(1,0,0).

    API kept compatible with prior usage:
      - length() -> number of meshes
      - generate_batch(mesh_idx, sampling_size=None) -> (data, label)
        * If sampling_size is provided, it's checked against meta.
    """

    def __init__(
        self,
        dataset_root: str,
        batch_size: int,
        meta_path: Optional[str] = None,
        patch_root: Optional[str] = None,
        drop_last: bool = True,
        shuffle_faces: bool = False,
        mmap: bool = True,
        rotation_anchor: str = 'center',  # 'center' or 'mean'
        eps: float = 1e-6,
    ):
        self.dataset_root = dataset_root
        self.patch_root = patch_root
        self.batch_size = int(batch_size)
        self.drop_last = bool(drop_last)
        self.shuffle_faces = bool(shuffle_faces)
        self.mmap_mode = 'r' if mmap else None
        self.rotation_anchor = rotation_anchor
        self.eps = eps

        # meta
        meta_path = meta_path or _default_meta_path(dataset_root)
        self.meta = _load_meta(meta_path)
        self.lsd_r_size = int(self.meta['lsd_r_size'])
        self.lsd_t_size = int(self.meta['lsd_t_size'])
        self.patch_num = int(self.meta['patch_num'])
        self.sampling_size = self.lsd_r_size * self.lsd_t_size + 1

        # mesh dirs
        self.mesh_dirs = _list_mesh_dirs(dataset_root)
        if not self.mesh_dirs:
            raise RuntimeError(f"No mesh dirs with lsd.npy & gt.npy under '{dataset_root}'")

        # Prepare a table of (lsd_path, gt_path, patch_path, nfaces)
        
        self._records: List[Tuple[str, str, str, int]] = []
        for mdir in self.mesh_dirs:
            lsd_path = os.path.join(mdir, 'lsd.npy')
            gt_path = os.path.join(mdir, 'gt.npy')
            patch_path = _resolve_patch_path(mdir, patch_root=self.patch_root, dataset_root=self.dataset_root)

            # Lazily open to validate shapes and collect nfaces
            lsd = np.load(lsd_path, mmap_mode=self.mmap_mode)
            gt = np.load(gt_path, mmap_mode=self.mmap_mode)
            patches = np.load(patch_path, mmap_mode=self.mmap_mode)

            # Validate shapes/dtypes
            if lsd.ndim != 3 or lsd.shape[1] != self.sampling_size or lsd.shape[2] != 3:
                raise ValueError(
                    f"{lsd_path}: expected (nfaces,{self.sampling_size},3), got {lsd.shape}")
            if gt.ndim != 2 or gt.shape[1] != 3:
                raise ValueError(f"{gt_path}: expected (nfaces,3), got {gt.shape}")
            if patches.ndim != 2 or patches.shape[0] != lsd.shape[0] or patches.shape[1] != self.patch_num:
                raise ValueError(
                    f"{patch_path}: expected (nfaces,{self.patch_num}), got {patches.shape}")

            nfaces = int(lsd.shape[0])
            self._records.append((lsd_path, gt_path, patch_path, nfaces))

    # ---- public API -----------------------------------------------------
    def length(self) -> int:
        return len(self._records)

    def get_meta_data(self) -> Dict[str, int]:
        return {
            'lsd_r_size': self.lsd_r_size,
            'lsd_t_size': self.lsd_t_size,
            'patch_num': self.patch_num,
        }

    def generate_batch(self, mesh_idx: int, sampling_size: Optional[int] = None):
        """
        Build per-mesh batches by iterating center faces -> patches.

        Returns:
            data:  (num_batches, B, M, N, 3)
            label: (num_batches, B, M, 3)
        """
        if sampling_size is not None and sampling_size != self.sampling_size:
            raise ValueError(
                f"sampling_size mismatch: got {sampling_size}, meta says {self.sampling_size}")

        lsd_path, gt_path, patch_path, nfaces = self._records[int(mesh_idx)]
        lsd = np.load(lsd_path, mmap_mode=self.mmap_mode)
        gt = np.load(gt_path, mmap_mode=self.mmap_mode)
        patch_faces = np.load(patch_path, mmap_mode=self.mmap_mode)

        # Sanity checks on indices
        if patch_faces.min() < 0 or patch_faces.max() >= nfaces:
            raise IndexError(
                f"Illegal face index in {patch_path}. Valid range [0,{nfaces-1}]")
        # Ensure center at first position (optional but handy)
        # If not guaranteed by generator, consider enforcing it here.

        # order of center faces
        center_idx = np.arange(nfaces, dtype=np.int64)
        if self.shuffle_faces:
            np.random.shuffle(center_idx)

        # batching over center faces
        B = self.batch_size
        total = len(center_idx)
        num_batches = total // B
        if not self.drop_last and total % B != 0:
            num_batches += 1

        batches_X: List[np.ndarray] = []
        batches_Y: List[np.ndarray] = []

        # target direction for rotation
        target = np.array([1.0, 0.0, 0.0], dtype=np.float64)

        for bi in range(num_batches):
            start = bi * B
            end = min(start + B, total)
            cur_centers = center_idx[start:end]
            cur_B = len(cur_centers)
            if cur_B < B and self.drop_last:
                break  # should not happen due to num_batches, but guard anyway

            # allocate batch arrays
            Xb = np.empty((cur_B, self.patch_num, self.sampling_size, 3), dtype=np.float32)
            Yb = np.empty((cur_B, self.patch_num, 3), dtype=np.float32)

            for i, c in enumerate(cur_centers):
                faces = patch_faces[c]  # (M,)
                # gather
                X = lsd[faces]  # (M,N,3)
                Y = gt[faces]   # (M,3)

                # patch-wise rotation
                if self.rotation_anchor == 'center':
                    anchor = Y[0]
                else:  # 'mean'
                    # mean of normalized Y to be robust
                    Yn = Y / np.maximum(np.linalg.norm(Y, axis=1, keepdims=True), self.eps)
                    anchor = Yn.mean(axis=0)

                if _safe_norm(anchor, self.eps) < 1e-6 or not np.isfinite(anchor).all():
                    # fallback: mean of non-degenerate gt
                    mask = np.isfinite(Y).all(axis=1) & (np.linalg.norm(Y, axis=1) > 1e-6)
                    if mask.any():
                        anchor = Y[mask].mean(axis=0)
                    else:
                        anchor = np.array([1.0, 0.0, 0.0], dtype=np.float64)

                R = _rotation_matrix_from_a_to_b(anchor, target, eps=self.eps)  # 3x3
                # apply same R to entire patch (X and Y)
                X_rot = X.reshape(-1, 3) @ R.T
                X_rot = X_rot.reshape(self.patch_num, self.sampling_size, 3)
                Y_rot = Y @ R.T

                # patch-level z-score on X
                Xn = _zscore_patch(X_rot, eps=self.eps)

                Xb[i] = Xn.astype(np.float32)
                Yb[i] = Y_rot.astype(np.float32)

            # if last batch smaller than B and drop_last=False, we keep cur_B
            batches_X.append(Xb)
            batches_Y.append(Yb)

        # stack to (num_batches, B, M, N, 3) and (num_batches, B, M, 3)
        if not batches_X:
            # no batches (e.g., nfaces < B and drop_last=True)
            return (
                np.empty((0, self.batch_size, self.patch_num, self.sampling_size, 3), dtype=np.float32),
                np.empty((0, self.batch_size, self.patch_num, 3), dtype=np.float32),
            )

        # Note: final batch may have size < B; caller should iterate first dim
        data = np.stack(batches_X, axis=0)
        label = np.stack(batches_Y, axis=0)
        return data, label
