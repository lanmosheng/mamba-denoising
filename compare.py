#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare THREE KPIs between two denoising result CSVs (from mesh_profile_eval.py).

KPI（越低越好）：
1) Macro mean of mean_deg (deg)
2) Count-weighted mean of mean_deg (deg)  —— 顶点模式等于 Micro；面模式是“网格大小加权”的近似
3) Macro mean of mse_mean (unit-normal vectors)

匹配规则（保证公平对比）：
  key = "{base}|{noise}|{gt_basename}"
  - base/noise 来自 pred 文件名（如 "armadillo_n1_01.off" -> base="armadillo", noise=1）
  - gt_basename 来自 GT 文件名（比如 "bunny_hi.obj"）

用法：
  python compare.py --new mamba_results.csv --base resnet_results.csv --labels "Mamba,ResNet"
"""
import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, help="CSV for NEW model (e.g., Mamba)")
    ap.add_argument("--base", required=True, help="CSV for BASELINE model (e.g., ResNet)")
    ap.add_argument("--labels", default="New,Baseline", help='Comma-separated labels, e.g. "Mamba,ResNet"')
    return ap.parse_args()


def _ok_rows(df: pd.DataFrame) -> pd.DataFrame:
    """保留成功的行：error 列为空/缺失视为成功；若无 error 列，则全部保留。"""
    if "error" not in df.columns:
        return df.copy()
    e = df["error"].astype(str)
    return df[(e.isna()) | (e.str.len() == 0) | (e.str.lower() == "nan") | (e.str.lower() == "none")].copy()


def parse_base_noise(pred_path: str) -> str:
    """
    从 pred 文件名解析 base 与 noise。
    例如: armadillo_n1_01.off -> base="armadillo", noise=1 -> 返回 "armadillo|1"
    若解析不到噪声等级，则退化为去掉扩展名的 basename。
    """
    name = Path(pred_path).name  # e.g., armadillo_n1_01.off
    m = re.match(r"(.+?)_n(\d)_[^\.]+", name)
    if m:
        base, noise = m.group(1), int(m.group(2))
        return f"{base}|{noise}"
    # fallback: use basename without extension
    return os.path.splitext(name)[0]


def make_key_row(row: pd.Series) -> str:
    gtname = Path(str(row.get("gt", ""))).name
    stem = parse_base_noise(str(row.get("pred", "")))
    return f"{stem}|{gtname}"


def prepare_df(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    ok = _ok_rows(df)
    ok["__key__"] = ok.apply(make_key_row, axis=1)
    return ok


def pick_cols(merged: pd.DataFrame, side: str, base_name: str, alt_name: str):
    """
    根据后缀选择列；优先 _new/_base，若无则兼容 _x/_y。
    side ∈ {"new","base"}；base_name/alt_name 是不带后缀的列名（如 "mean_deg"）
    """
    if side == "new":
        prefer = f"{base_name}_new"
        alt = f"{base_name}_x"
    else:
        prefer = f"{base_name}_base"
        alt = f"{base_name}_y"

    if prefer in merged.columns:
        return prefer
    if alt in merged.columns:
        return alt
    # 若两者都不存在（例如旧CSV缺少 mse_mean），返回 None
    return None


def kpis_from_side(merged: pd.DataFrame, side: str):
    """
    返回 (macro_angle_mean, count_weighted_angle_mean, macro_mse_mean)
    若 CSV 无 mse_mean，将返回 macro_mse_mean = nan，并在外层提示。
    """
    mean_col = pick_cols(merged, side, "mean_deg", "mean_deg")
    cnt_col  = pick_cols(merged, side, "count", "count")
    mse_col  = pick_cols(merged, side, "mse_mean", "mse_mean")  # 可能不存在（老CSV）

    if mean_col is None or cnt_col is None:
        raise KeyError(f"Cannot find expected columns for side='{side}'. Have: {list(merged.columns)}")

    angle_vals = merged[mean_col].to_numpy(dtype=float)
    counts = merged[cnt_col].to_numpy(dtype=float)

    macro_angle = float(np.nanmean(angle_vals))
    cw_angle = float(np.average(angle_vals, weights=counts)) if np.isfinite(counts).all() and counts.sum() > 0 else float("nan")

    macro_mse = float("nan")
    if mse_col is not None:
        mse_vals = merged[mse_col].to_numpy(dtype=float)
        macro_mse = float(np.nanmean(mse_vals))

    return macro_angle, cw_angle, macro_mse


def fmt(x) -> str:
    return "N/A" if (x is None or (isinstance(x, float) and not np.isfinite(x))) else f"{x:.4f}"


def rel_improv(new: float, base: float) -> float:
    """相对提升百分比： (base - new) / base * 100。base 为 0 或无效返回 NaN。"""
    if base is None or not np.isfinite(base) or base == 0:
        return float("nan")
    if new is None or not np.isfinite(new):
        return float("nan")
    return (base - new) / base * 100.0


def main():
    args = parse_args()
    label_new, label_base = [s.strip() for s in args.labels.split(",", 1)]

    a = prepare_df(args.new)
    b = prepare_df(args.base)

    # 合并时强制指定后缀，避免 _x/_y 的混乱
    merged = a.merge(b, on="__key__", suffixes=("_new", "_base"))
    matched = len(merged)
    total_new = len(a)
    total_base = len(b)

    if matched == 0:
        # 提示可能的 key 差异
        keys_a = set(a["__key__"].tolist())
        keys_b = set(b["__key__"].tolist())
        only_a = list(keys_a - keys_b)[:5]
        only_b = list(keys_b - keys_a)[:5]
        raise SystemExit(
            "No matched pairs between CSVs.\n"
            f"- Example keys only in NEW ({len(keys_a-keys_b)}): {only_a}\n"
            f"- Example keys only in BASE ({len(keys_b-keys_a)}): {only_b}\n"
            "Check file naming pattern or GT filenames."
        )

    # 计算三项 KPI（各自）
    m_mean_new, cw_new, m_mse_new = kpis_from_side(merged, "new")
    m_mean_base, cw_base, m_mse_base = kpis_from_side(merged, "base")

    # 计算差与相对提升（正数=新模型更好）
    delta_m = m_mean_new - m_mean_base
    delta_cw = cw_new - cw_base
    delta_mse = m_mse_new - m_mse_base

    ri_m = rel_improv(m_mean_new, m_mean_base)
    ri_cw = rel_improv(cw_new, cw_base)
    ri_mse = rel_improv(m_mse_new, m_mse_base)

    print("=== 3-KPI Comparison ===")
    print(f"Matched pairs: {matched} (new {total_new}, base {total_base})\n")

    print("1) Macro mean of mean_deg (deg)")
    print(f"   {label_new}: {fmt(m_mean_new)}   | {label_base}: {fmt(m_mean_base)}   "
          f"| Δ(new-base): {fmt(delta_m)}   | Rel. improv: {fmt(ri_m)}%")

    print("2) Count-weighted mean of mean_deg (deg)")
    print(f"   {label_new}: {fmt(cw_new)}   | {label_base}: {fmt(cw_base)}   "
          f"| Δ(new-base): {fmt(delta_cw)}   | Rel. improv: {fmt(ri_cw)}%")

    print("3) Macro mean of MSE (unit-normal vectors)")
    print(f"   {label_new}: {fmt(m_mse_new)}   | {label_base}: {fmt(m_mse_base)}   "
          f"| Δ(new-base): {fmt(delta_mse)}   | Rel. improv: {fmt(ri_mse)}%")

    if not np.isfinite(m_mse_new) or not np.isfinite(m_mse_base):
        print("\n[Note] CSV 中未找到 mse_mean 列，第三项 KPI 显示为 N/A。"
              " 请使用更新后的 mesh_profile_eval.py 生成带 MSE 的 CSV。")


if __name__ == "__main__":
    main()
