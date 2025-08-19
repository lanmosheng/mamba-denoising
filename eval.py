#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch mesh denoising evaluation with ANGLE and MSE metrics.

PROFILE FORMAT
--------------
Line 1: an integer N = number of pairs
Then repeat N times (exactly 2 lines per pair):
  line: <pred_mesh_path>
  line: <gt_mesh_path>

Blank lines and lines starting with '#' are ignored.

USAGE
-----
# 与训练一致（方向敏感，非 acute），面模式 + 面积加权，输出 CSV
python mesh_profile_eval.py --profile eval.profile --mode face --weighted --csv results.csv

# 若需要旧报告口径（方向不敏感），加 --acute
python mesh_profile_eval.py --profile eval.profile --mode face --weighted --acute --csv results.csv

DEPENDENCIES
------------
pip install numpy trimesh
"""
import argparse
import csv
import os
import sys
from typing import List, Tuple, Dict, Optional

import numpy as np
import trimesh


def clamp(v, vmin=-1.0, vmax=1.0):
    return np.minimum(np.maximum(v, vmin), vmax)


def angle_deg_from_dot(dot):
    dot = clamp(dot)
    return np.degrees(np.arccos(dot))


def face_normals_and_areas(mesh: trimesh.Trimesh):
    fn = mesh.face_normals
    areas = mesh.area_faces
    return fn, areas


def vertex_normals(mesh: trimesh.Trimesh):
    return mesh.vertex_normals


def load_single_mesh(path: str) -> trimesh.Trimesh:
    obj = trimesh.load(path, process=True)
    if isinstance(obj, trimesh.Scene):
        meshes = [g for g in obj.dump() if isinstance(g, trimesh.Trimesh)]
        if not meshes:
            raise ValueError(f"No triangle geometry in scene: {path}")
        mesh = trimesh.util.concatenate(meshes)
    elif isinstance(obj, trimesh.Trimesh):
        mesh = obj
    else:
        raise ValueError(f"Unsupported geometry type for: {path}")
    return mesh


def stats(angles_deg: np.ndarray, weights: Optional[np.ndarray] = None):
    if weights is not None:
        w = weights / (np.sum(weights) + 1e-12)
        mean = float(np.sum(angles_deg * w))
    else:
        mean = float(np.mean(angles_deg))
    out = {
        "mean_deg": mean,
        "median_deg": float(np.median(angles_deg)),
        "std_deg": float(np.std(angles_deg)),
        "p90_deg": float(np.percentile(angles_deg, 90)),
        "p95_deg": float(np.percentile(angles_deg, 95)),
        "max_deg": float(np.max(angles_deg)),
        "count": int(angles_deg.size),
    }
    return out


def stats_mse(mse_vals: np.ndarray, weights: Optional[np.ndarray] = None):
    if weights is not None:
        w = weights / (np.sum(weights) + 1e-12)
        mean = float(np.sum(mse_vals * w))
    else:
        mean = float(np.mean(mse_vals))
    out = {
        "mse_mean": mean,
        "mse_median": float(np.median(mse_vals)),
        "mse_std": float(np.std(mse_vals)),
        "mse_p90": float(np.percentile(mse_vals, 90)),
        "mse_p95": float(np.percentile(mse_vals, 95)),
        "mse_max": float(np.max(mse_vals)),
    }
    return out


def compute_pair(pred_path: str, gt_path: str, mode: str, acute: bool, weighted: bool
                 ) -> Tuple[Dict[str, float], float, float, float, int]:
    """
    Returns:
      stats_dict (angles + mse),
      sum_aw_angle, sum_w, sum_aw_mse, count_samples
    """
    pred = load_single_mesh(pred_path)
    gt = load_single_mesh(gt_path)

    if mode == "face":
        if pred.faces.shape[0] != gt.faces.shape[0]:
            raise ValueError(f"Face counts differ: pred={pred.faces.shape[0]} vs gt={gt.faces.shape[0]}")
        fn_pred, area_pred = face_normals_and_areas(pred)
        fn_gt, _ = face_normals_and_areas(gt)
        # Unit normalization (safety)
        fn_pred = fn_pred / (np.linalg.norm(fn_pred, axis=1, keepdims=True) + 1e-12)
        fn_gt = fn_gt / (np.linalg.norm(fn_gt, axis=1, keepdims=True) + 1e-12)
        dots = np.sum(fn_pred * fn_gt, axis=1)
        if acute:
            dots = np.abs(dots)
        ang = angle_deg_from_dot(dots)                         # (F,)
        mse = np.sum((fn_pred - fn_gt) ** 2, axis=1)           # (F,)
        w = area_pred.astype(np.float64) if weighted else np.ones_like(ang, dtype=np.float64)
        use_weights = w if weighted else None                  # 仅 face+--weighted 下对 mean 加权
    else:  # vertex
        if pred.vertices.shape[0] != gt.vertices.shape[0]:
            raise ValueError(f"Vertex counts differ: pred={pred.vertices.shape[0]} vs gt={gt.vertices.shape[0]}")
        vn_pred = vertex_normals(pred)
        vn_gt = vertex_normals(gt)
        vn_pred = vn_pred / (np.linalg.norm(vn_pred, axis=1, keepdims=True) + 1e-12)
        vn_gt = vn_gt / (np.linalg.norm(vn_gt, axis=1, keepdims=True) + 1e-12)
        dots = np.sum(vn_pred * vn_gt, axis=1)
        if acute:
            dots = np.abs(dots)
        ang = angle_deg_from_dot(dots)                         # (V,)
        mse = np.sum((vn_pred - vn_gt) ** 2, axis=1)           # (V,)
        w = np.ones_like(ang, dtype=np.float64)
        use_weights = None                                     # vertex 模式不加权 mean

    sum_w = float(np.sum(w))
    sum_aw_angle = float(np.sum(ang * w))
    sum_aw_mse = float(np.sum(mse * w))

    angle_stats = stats(ang, use_weights)
    mse_stats = stats_mse(mse, use_weights)

    row_stats = {**angle_stats, **mse_stats}
    return row_stats, sum_aw_angle, sum_w, sum_aw_mse, int(ang.size)


def parse_profile(path: str):
    with open(path, 'r', encoding='utf-8') as f:
        raw_lines = [ln.strip() for ln in f.readlines()]
    lines = [ln for ln in raw_lines if ln and not ln.lstrip().startswith('#')]
    if not lines:
        raise ValueError("Empty profile.")
    try:
        n = int(lines[0])
    except Exception:
        raise ValueError("First non-comment line must be an integer N.")
    pairs = []
    idx = 1
    for i in range(n):
        if idx + 1 >= len(lines):
            raise ValueError(f"Insufficient file paths for pair {i}. Expected 2 lines after count.")
        pred = lines[idx]; gt = lines[idx + 1]
        pairs.append((pred, gt))
        idx += 2
    return pairs


def main():
    ap = argparse.ArgumentParser(description="Batch evaluate mesh angle error and MSE using a profile listing mesh pairs.")
    ap.add_argument("--profile", required=True, help="Path to profile file.")
    ap.add_argument("--mode", choices=["face", "vertex"], default="face",
                    help="Compare per-face or per-vertex normals. Default: face")
    ap.add_argument("--acute", action="store_true",
                    help="Use acute angle acos(|dot|). Omit for training-consistent, direction-sensitive evaluation.")
    ap.add_argument("--weighted", action="store_true",
                    help="Area-weighted mean (face mode only). Ignored in vertex mode.")
    ap.add_argument("--csv", default="",
                    help="Optional path to write CSV summary. Default: <profile>_angle_results.csv")
    ap.add_argument("--strict", action="store_true",
                    help="Stop on first error (default: continue and record error).")
    args = ap.parse_args()

    try:
        pairs = parse_profile(args.profile)
    except Exception as e:
        print(f"[Error] Bad profile: {e}", file=sys.stderr)
        sys.exit(2)

    results = []

    # Global sums for micro-averages
    tot_w = 0.0
    tot_aw_angle = 0.0
    tot_aw_mse = 0.0

    num_ok = 0
    num_err = 0

    for i, (pred, gt) in enumerate(pairs):
        row = {
            "index": i,
            "pred": pred,
            "gt": gt,
            "mode": args.mode,
            "acute": bool(args.acute),
            "weighted": bool(args.weighted and args.mode == "face"),
        }
        try:
            stats_dict, sum_aw_ang, sum_w, sum_aw_mse, cnt = compute_pair(
                pred, gt, args.mode, args.acute, args.weighted
            )
            row.update(stats_dict)
            row["error"] = ""
            # Accumulate micro sums
            tot_w += sum_w
            tot_aw_angle += sum_aw_ang
            tot_aw_mse += sum_aw_mse
            num_ok += 1
        except Exception as e:
            row.update({
                "mean_deg": np.nan, "median_deg": np.nan, "std_deg": np.nan,
                "p90_deg": np.nan, "p95_deg": np.nan, "max_deg": np.nan, "count": 0,
                "mse_mean": np.nan, "mse_median": np.nan, "mse_std": np.nan,
                "mse_p90": np.nan, "mse_p95": np.nan, "mse_max": np.nan,
                "error": str(e),
            })
            num_err += 1
            if args.strict:
                print(f"[Error] Pair {i} failed: {e}", file=sys.stderr)
                print("Stopping due to --strict.")
                results.append(row)
                break
        results.append(row)

    print("=== Batch Evaluation Summary (training-consistent by default: NON-ACUTE) ===")
    print(f"Pairs total: {len(pairs)} | OK: {num_ok} | Error: {num_err}")
    if tot_w > 0:
        micro_mean_angle = tot_aw_angle / tot_w
        micro_mean_mse = tot_aw_mse / tot_w
        print(f"Micro-averaged ANGLE mean (deg): {micro_mean_angle:.6f}")
        print(f"Micro-averaged MSE (unit-normal vectors): {micro_mean_mse:.6f}")

    if num_ok > 0:
        macro_mean_angle = np.nanmean([r["mean_deg"] for r in results if r["error"] == ""])
        macro_mean_mse = np.nanmean([r["mse_mean"] for r in results if r["error"] == ""])
        print(f"Macro-averaged ANGLE mean (deg): {macro_mean_angle:.6f}")
        print(f"Macro-averaged MSE (unit-normal vectors): {macro_mean_mse:.6f}")

    # CSV output
    if args.csv:
        outp = args.csv
    else:
        base = os.path.splitext(os.path.basename(args.profile))[0]
        outp = f"{base}_angle_results.csv"

    fieldnames = ["index", "pred", "gt", "mode", "acute", "weighted",
                  # angle stats
                  "mean_deg", "median_deg", "std_deg", "p90_deg", "p95_deg", "max_deg", "count",
                  # mse stats
                  "mse_mean", "mse_median", "mse_std", "mse_p90", "mse_p95", "mse_max",
                  "error"]
    try:
        with open(outp, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow({k: r.get(k, "") for k in fieldnames})
        print(f"Wrote CSV: {outp}")
    except Exception as e:
        print(f"[Warning] Failed to save CSV: {e}", file=sys.stderr)

    # Compact per-row summary
    for r in results:
        if r["error"]:
            print(f"[{r['index']}] ERROR: {r['error']} | pred={r['pred']} | gt={r['gt']}")
        else:
            print(f"[{r['index']}] mean={r['mean_deg']:.6f}° | p90={r['p90_deg']:.3f}° | max={r['max_deg']:.3f}°"
                  f" | mse_mean={r['mse_mean']:.6f} | cnt={r['count']}"
                  f" | pred={r['pred']} | gt={r['gt']}")
    

if __name__ == "__main__":
    main()
