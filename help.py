#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import json

# 路径改成你要检查的 LSD 根目录
LSD_ROOT = "valset/data1001"

def check_lsd_shapes(lsd_root):
    meta_path = os.path.join(lsd_root, "meta.json")
    if not os.path.isfile(meta_path):
        raise FileNotFoundError(f"meta.json not found: {meta_path}")

    # 读 meta.json，得到标准采样大小
    with open(meta_path, "r", encoding="utf-8") as f:
        meta = json.load(f)
    r, t = meta["lsd_r_size"], meta["lsd_t_size"]
    sampling_size = r * t + 1

    print(f"[info] Expect shape = (nF, {sampling_size}, 3)")

    bad_files = []
    mesh_dirs = [d for d in os.listdir(lsd_root) if os.path.isdir(os.path.join(lsd_root, d))]

    for name in sorted(mesh_dirs):
        lsd_path = os.path.join(lsd_root, name, "lsd.npy")
        if not os.path.isfile(lsd_path):
            print(f"[skip] {name}: lsd.npy not found")
            continue
        try:
            arr = np.load(lsd_path)
        except Exception as e:
            print(f"[error] {name}: cannot load ({e})")
            bad_files.append((name, "load_error"))
            continue

        if arr.ndim != 3 or arr.shape[1] != sampling_size or arr.shape[2] != 3:
            print(f"[bad] {name}: shape={arr.shape}")
            bad_files.append((name, arr.shape))
        else:
            print(f"[ok] {name}: shape={arr.shape}")

    print("\n===== Summary =====")
    if bad_files:
        print(f"Found {len(bad_files)} bad LSD files:")
        for name, shape in bad_files:
            print(f"  - {name}: shape={shape}")
    else:
        print("All LSD files are valid.")

if __name__ == "__main__":
    check_lsd_shapes(LSD_ROOT)
