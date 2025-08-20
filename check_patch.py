# check_patch.py
import sys
import numpy as np
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_patch.py <patch_dir> [lsd_dir_or_nfaces] [expected_M]")
        sys.exit(2)

    patch_dir = Path(sys.argv[1])
    patch_path = patch_dir / "patch_faces.npy"
    if not patch_path.exists():
        print(f"❌ missing {patch_path}")
        sys.exit(1)

    patch = np.load(patch_path, mmap_mode="r")   # (nfaces, M)
    print(f"patch: shape={patch.shape}, dtype={patch.dtype}")
    nfaces, M = patch.shape

    # 读取第二个参数：可传 lsd 目录（含 lsd.npy），或直接传 nfaces 数字
    nfaces_expected = None
    if len(sys.argv) >= 3:
        p = Path(sys.argv[2])
        if p.is_dir() and (p / "lsd.npy").exists():
            lsd = np.load(p / "lsd.npy", mmap_mode="r")
            nfaces_expected = lsd.shape[0]
            print(f"lsd: shape={lsd.shape}, dtype={lsd.dtype}")
        else:
            nfaces_expected = int(sys.argv[2])

    expected_M = int(sys.argv[3]) if len(sys.argv) >= 4 else None

    # 形状与类型
    if nfaces_expected is not None and nfaces != nfaces_expected:
        print(f"❌ nfaces mismatch: patch={nfaces} vs expected={nfaces_expected}")
        sys.exit(1)
    if expected_M is not None and M != expected_M:
        print(f"❌ M mismatch: patch={M} vs expected={expected_M}")
        sys.exit(1)
    if patch.dtype != np.int32:
        print(f"⚠️ dtype is {patch.dtype}, recommended int32")

    # 索引范围检查
    bad = (patch < 0) | (patch >= nfaces)
    if bad.any():
        loc = np.argwhere(bad)
        print(f"❌ out-of-range indices (showing first 10): {loc[:10].tolist()} "
              f"total={bad.sum()}")
        sys.exit(1)
    else:
        print("✅ all indices in [0, nfaces)")

    # 第一列是否中心面（如果你约定 center_is_first）
    center_rate = (patch[:, 0] == np.arange(nfaces)).mean()
    print(f"center_is_first rate: {center_rate*100:.2f}%")

    # 估计 padding 情况：统计前 1w 行的唯一面数
    sample = min(nfaces, 10_000)
    uniq_counts = np.array([len(np.unique(row)) for row in patch[:sample]])
    padded_ratio = (uniq_counts < M).mean() * 100
    print(f"rows with padding (first {sample}): {padded_ratio:.2f}%")

    # 覆盖度：每个面被包含次数
    freq = np.bincount(patch.ravel(), minlength=nfaces)
    print(f"coverage per-face: min={freq.min()}, mean={freq.mean():.1f}, max={freq.max()}")
    never = np.where(freq == 0)[0]
    if len(never):
        print(f"❌ faces never included in any patch: count={len(never)} "
              f"examples={never[:10].tolist()}")
    else:
        print("✅ every face appears in at least one patch")

    print("Done.")

if __name__ == "__main__":
    main()
