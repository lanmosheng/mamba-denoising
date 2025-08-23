# -*- coding: utf-8 -*-
# 从 patch_root/{mesh_name}/patch_faces.npy 读取，做 IOU 去重，输出 centers.npy
# 终端打印：usable_patches / coverage / avg_cov / max_cov / coverage_hist
import argparse, random
from pathlib import Path
import numpy as np

# 尝试用 tqdm；没装则优雅降级为无进度条
try:
    from tqdm import tqdm
except Exception:
    class tqdm:
        def __init__(self, iterable=None, total=None, desc=None, leave=False):
            self.iterable = iterable if iterable is not None else range(total or 0)
        def __iter__(self):
            for x in self.iterable:
                yield x
        def set_postfix_str(self, s): pass
        def close(self): pass

def list_meshes_from_patch_root(patch_root: Path):
    return [d for d in sorted(patch_root.iterdir()) if d.is_dir() and (d / "patch_faces.npy").exists()]

def jaccard_iou(a: np.ndarray, b: np.ndarray) -> float:
    sa, sb = set(map(int, a)), set(map(int, b))
    inter = len(sa & sb)
    union = len(sa) + len(sb) - inter
    return inter / union if union else 0.0

def build_centers_for_mesh(patch_faces: np.ndarray, iou_thr: float,
                           shuffle: bool, seed: int, max_centers: int | None,
                           use_pbar: bool, mesh_name: str) -> np.ndarray:
    nfaces = patch_faces.shape[0]
    order = list(range(nfaces))
    if shuffle:
        rnd = random.Random(seed)
        rnd.shuffle(order)

    iterator = tqdm(order, desc=f"select[{mesh_name}]", leave=False) if use_pbar else order
    selected: list[int] = []
    for idx in iterator:
        p = patch_faces[idx]
        keep = True
        # 与已选中心比 IOU（贪心去重）
        for s in selected:
            if jaccard_iou(p, patch_faces[s]) >= iou_thr:
                keep = False
                break
        if keep:
            selected.append(idx)
            if use_pbar and hasattr(iterator, "set_postfix_str"):
                iterator.set_postfix_str(f"K={len(selected)}")
            if max_centers is not None and len(selected) >= max_centers:
                break
    if use_pbar and hasattr(iterator, "close"):
        iterator.close()
    return np.array(selected, dtype=np.int64)

def coverage_stats(patch_faces: np.ndarray, centers: np.ndarray):
    nfaces = patch_faces.shape[0]
    cov = np.zeros(nfaces, dtype=np.int32)
    for c in centers:
        cov[patch_faces[c]] += 1
    covered = int((cov > 0).sum())
    coverage_rate = covered / nfaces if nfaces else 0.0
    avg_cov = float(cov.mean()) if nfaces else 0.0
    max_cov = int(cov.max()) if cov.size else 0
    hist = np.bincount(cov) if cov.size else np.array([])
    hist_str = " ".join(f"{k}:{int(v)}" for k, v in enumerate(hist[:30])) if hist.size else ""
    return coverage_rate, covered, avg_cov, max_cov, hist_str

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--patch-root', required=True, help='例如 patches/')
    ap.add_argument('--iou', type=float, default=0.6)
    ap.add_argument('--shuffle', action='store_true', help='先随机打乱候选面顺序再做贪心')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--max-centers', type=int, default=None)
    ap.add_argument('--no-pbar', action='store_true', help='关闭进度条（更快）')
    args = ap.parse_args()

    patch_root = Path(args.patch_root)
    mesh_dirs = list_meshes_from_patch_root(patch_root)
    if not mesh_dirs:
        raise SystemExit(f"No patch_faces.npy found under {patch_root}")

    outer = tqdm(mesh_dirs, desc="meshes", leave=False) if not args.no_pbar else mesh_dirs
    for mdir in outer:
        mesh_name = mdir.name
        patch_path = mdir / 'patch_faces.npy'
        patches = np.load(patch_path, mmap_mode='r')  # (nfaces, M)
        nfaces, M = patches.shape

        centers = build_centers_for_mesh(
            patches, args.iou, args.shuffle, args.seed, args.max_centers,
            use_pbar=(not args.no_pbar), mesh_name=mesh_name
        )

        cov_rate, covered, avg_cov, max_cov, hist_str = coverage_stats(patches, centers)
        usable = len(centers)

        # 保存 centers.npy
        np.save(mdir / 'centers.npy', centers)

        # 终端输出
        print(
            f"[OK] {mesh_name} | nfaces={nfaces} | usable_patches={usable} "
            f"| coverage={cov_rate:.3f} ({covered}/{nfaces}) "
            f"| avg_cov={avg_cov:.2f} | max_cov={max_cov}"
        )
        print(f"     coverage_hist: {hist_str}")

if __name__ == '__main__':
    main()
