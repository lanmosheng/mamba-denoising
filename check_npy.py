# check_x_bias.py
import sys, numpy as np
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_x_bias.py <mesh_dir> [threshold_deg]")
        sys.exit(2)

    mesh_dir = Path(sys.argv[1])
    th_deg = float(sys.argv[2]) if len(sys.argv) >= 3 else 11.5  # 约等于 cos>0.98
    cos_th = np.cos(np.deg2rad(th_deg))

    lsd = np.load(mesh_dir/"lsd.npy", mmap_mode="r")   # (nfaces, N, 3)
    gt  = np.load(mesh_dir/"gt.npy",  mmap_mode="r")   # (nfaces, 3)

    ex = np.array([1.0, 0.0, 0.0], dtype=np.float64)

    # --- GT 偏向 ---
    gt64 = gt.astype(np.float64)
    gt_norm = np.linalg.norm(gt64, axis=1, keepdims=True)
    gt_norm[gt_norm==0] = 1.0
    gt_cos = (gt64/gt_norm) @ ex
    frac_gt_aligned = float((gt_cos > cos_th).mean())
    print(f"GT aligned to +X (> {th_deg}°): {frac_gt_aligned*100:.2f}%")

    # --- LSD 偏向（全部 token）---
    nfaces, N, _ = lsd.shape
    A = lsd.reshape(-1, 3).astype(np.float64)
    A_norm = np.linalg.norm(A, axis=1, keepdims=True)
    A_norm[A_norm==0] = 1.0
    lsd_cos = (A/A_norm) @ ex
    frac_lsd_aligned = float((lsd_cos > cos_th).mean())
    print(f"LSD aligned to +X (> {th_deg}°): {frac_lsd_aligned*100:.2f}%")

    # --- LSD 每个面的“近中心层”粗估（可选）：取每面前 K 个 token 看偏向 ---
    # 若你的采样顺序是按半径从小到大，设置一个小 K（如 1~3）更接近“中心半径”。
    K = min(3, N)
    early = lsd[:, :K, :].reshape(-1, 3).astype(np.float64)
    e_norm = np.linalg.norm(early, axis=1, keepdims=True); e_norm[e_norm==0] = 1.0
    frac_lsd_early = float(((early/e_norm) @ ex > cos_th).mean())
    print(f"LSD first {K} tokens aligned to +X: {frac_lsd_early*100:.2f}%")

    # 粗判定：如果这些占比非常高（比如 >50%），基本就是“还在对齐到 +X”
    warn = []
    if frac_gt_aligned > 0.5: warn.append("GT strong +X bias")
    if frac_lsd_aligned > 0.5: warn.append("LSD strong +X bias")
    if frac_lsd_early > 0.5: warn.append("LSD-center strong +X bias")

    if warn:
        print("❌ Likely still rotating per-face:", "; ".join(warn))
    else:
        print("✅ No obvious +X bias (looks global-frame).")

if __name__ == "__main__":
    main()
