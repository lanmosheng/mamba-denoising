# generate_stage1_offline.py
# -*- coding: utf-8 -*-
import os, json, math, time
from typing import List, Dict, Tuple, Optional
import numpy as np
import torch
import torch.nn.functional as F
from tqdm import tqdm

# =========================
# 路径 & 生成配置（按需修改）
# =========================
DATASET_ROOT = "dataset_s_i1/train1601"     # mesh 根目录（包含若干 <mesh_name>/lsd.npy, gt.npy）
PATCH_ROOT   = "patches"                    # patch_faces.npy 根目录：<PATCH_ROOT>/<mesh_name>/patch_faces.npy
OUT_ROOT     = "offline_s1/train1601"       # 离线输出根目录（每个 mesh 一个子目录）

BEST_MODEL_PATH = "out/s1_mamba_face_1601/stage1_best.pt"  # 你训练出的 S1 best ckpt

# 模型构造超参（需与训练时一致）
D_MODEL     = 64
FACE_DEPTH  = 4
D_STATE     = 16
D_CONV      = 4
EXPAND      = 2
RES_SCALE   = 0.5
P_DROP      = 0.05

# 执行 & 性能配置
DEVICE                = torch.device("cuda" if torch.cuda.is_available() else "cpu")
USE_DATAPARALLEL      = (torch.cuda.device_count() >= 2)   # 多卡推理
ROW_CHUNK             = 64     # 每次处理多少行 patch（旋转/准备）
MICRO_BATCH_FACES     = 240   # 每次前向多少个“面”（R*M 会切块）
SAVE_N1               = True   # 是否同时保存 S1 的法向初稿 N1_patch
DTYPE_FEAT            = np.float16
DTYPE_Y               = np.float16  # 可改 np.float32

# =========================
# 工具函数（与训练/loader一致）
# =========================
def _list_mesh_dirs(dataset_root: str) -> List[str]:
    mesh_dirs = []
    for root, _, files in os.walk(dataset_root):
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
    if patch_root is None:
        raise ValueError("patch_root is required")
    mesh_name = os.path.basename(mesh_dir.rstrip(os.sep))
    cand = os.path.join(patch_root, mesh_name, 'patch_faces.npy')
    if os.path.exists(cand):
        return cand
    raise FileNotFoundError(f"patch_faces.npy not found for mesh '{mesh_name}' under patch_root='{patch_root}'")

def _rotation_matrix_a2b(a: np.ndarray, b: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """把单位向量 a 旋到 b 的 3x3 旋转矩阵（稳定处理平行/反平行）。"""
    v = np.cross(a, b)
    c = float(np.dot(a, b))    # cosθ
    s = np.linalg.norm(v)      # |v| = sinθ
    if s < eps:
        if c > 0.0:
            return np.eye(3, dtype=np.float64)  # a≈b
        # a≈-b，取任一不共线轴做 180°
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

CANONICAL_DIR = np.array([1.0, 0.0, 0.0], dtype=np.float64)

def _rotation_from_center_token0(X_patch: np.ndarray) -> np.ndarray:
    """
    用中心面（faces[0]）的 LSD token0 作为朝向向量 a，计算 R 使 a -> (1,0,0)。
    仅依赖 LSD。
    """
    a = X_patch[0, 0].astype(np.float64)
    a_norm = np.linalg.norm(a)
    if a_norm < 1e-12:
        return np.eye(3, dtype=np.float64)
    a /= a_norm
    return _rotation_matrix_a2b(a, CANONICAL_DIR)

# =========================
# 模型构建（与你项目保持一致）
# =========================
from train_utils import config  # 需要你项目里提供 get_face_mamba

def build_face_model(meta: Dict):
    r = int(meta['lsd_r_size'])
    t = int(meta['lsd_t_size'])
    model = config.get_face_mamba(
        DEVICE,
        lsd_r_size=r,
        lsd_t_size=t,
        d_model=D_MODEL,
        depth=FACE_DEPTH,
        d_state=D_STATE,
        d_conv=D_CONV,
        expand=EXPAND,
        residual_scale=RES_SCALE,
        p_drop=P_DROP,
    )
    # 加载 best 权重
    ckpt = torch.load(BEST_MODEL_PATH, map_location='cpu')
    (model.module if isinstance(model, torch.nn.DataParallel) else model).load_state_dict(ckpt['model'], strict=True)

    if USE_DATAPARALLEL and not isinstance(model, torch.nn.DataParallel) and DEVICE.type == "cuda":
        model = torch.nn.DataParallel(model, device_ids=list(range(torch.cuda.device_count())))
    model.eval()
    return model

# =========================
# 主流程：对每个 mesh 生成离线产物
# =========================
def process_one_mesh(mesh_dir: str, patch_root: str, out_dir: str,
                     model: torch.nn.Module, meta: Dict):
    os.makedirs(out_dir, exist_ok=True)

    lsd_path = os.path.join(mesh_dir, "lsd.npy")
    gt_path  = os.path.join(mesh_dir, "gt.npy")
    pf_path  = _resolve_patch_path(mesh_dir, PATCH_ROOT)

    lsd = np.load(lsd_path, mmap_mode='r')             # (nfaces, N, 3)
    gt  = np.load(gt_path,  mmap_mode='r')             # (nfaces, 3)
    pf  = np.load(pf_path,  mmap_mode='r')             # (K, M)
    K, M = pf.shape
    N    = 1 + int(meta['lsd_r_size']) * int(meta['lsd_t_size'])

    # 校验
    if lsd.shape[1] != N or lsd.shape[2] != 3:
        raise ValueError(f"{lsd_path}: expect (*,{N},3), got {lsd.shape}")
    if gt.shape[1] != 3:
        raise ValueError(f"{gt_path}: expect (*,3), got {gt.shape}")
    if pf.min() < 0 or pf.max() >= lsd.shape[0]:
        raise IndexError(f"{pf_path}: face index out of range")
    if int(meta['patch_num']) != M:
        print(f"[warn] meta.patch_num={meta['patch_num']} but patch_faces has M={M}")

    # 目标 memmap
    feat_path = os.path.join(out_dir, "F_patch.fp16.memmap")
    yrot_path = os.path.join(out_dir, "Y_patch_rot.fp16.memmap")
    n1_path   = os.path.join(out_dir, "N1_patch.fp16.memmap")

    Fmm = np.memmap(feat_path, mode='w+', dtype=DTYPE_FEAT, shape=(K, M, D_MODEL))
    Ymm = np.memmap(yrot_path, mode='w+', dtype=DTYPE_Y,   shape=(K, M, 3))
    Nmm = None
    if SAVE_N1:
        Nmm = np.memmap(n1_path, mode='w+', dtype=np.float16, shape=(K, M, 3))

    # 逐块处理
    total_rows = K
    pbar = tqdm(total=total_rows, desc=os.path.basename(mesh_dir), ncols=60, bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}")
    with torch.no_grad():
        for beg in range(0, K, ROW_CHUNK):
            end = min(K, beg + ROW_CHUNK)
            rows = pf[beg:end]                      # [R, M]
            R_list: List[np.ndarray] = []
            X_list: List[np.ndarray] = []
            Y_list: List[np.ndarray] = []

            # 先旋转到统一坐标系（与训练一致）
            for faces in rows:
                X_patch = lsd[faces]                # (M, N, 3)
                R = _rotation_from_center_token0(X_patch)  # (3,3)
                X_rot = X_patch @ R.T               # (M, N, 3)
                Y_rot = gt[faces] @ R.T             # (M, 3)
                R_list.append(R)
                X_list.append(X_rot.astype(np.float32))
                Y_list.append(Y_rot.astype(np.float32))

            # 合并成 (R*M, N, 3)，做 S1 前向
            XRM = np.concatenate(X_list, axis=0)    # (R*M, N, 3)
            X_tensor = torch.from_numpy(XRM).to(DEVICE, non_blocking=True)

            # 切成若干 micro-batch（避免显存爆）
            feats_out = np.empty((XRM.shape[0], D_MODEL), dtype=np.float32)
            n1_out    = np.empty((XRM.shape[0], 3),       dtype=np.float32) if SAVE_N1 else None

            for mb_beg in range(0, XRM.shape[0], MICRO_BATCH_FACES):
                mb_end = min(XRM.shape[0], mb_beg + MICRO_BATCH_FACES)
                xb = X_tensor[mb_beg:mb_end]       # (B, N, 3)
                out = model(xb)
                if isinstance(out, tuple):
                    n_pred, feat = out
                else:
                    # 你的 FaceEncoder 应该返回 (n, feat)。若不是，抛错更安全。
                    raise RuntimeError("FaceEncoder.forward 需要返回 (n, feat) 以便离线导出特征。")
                feats_out[mb_beg:mb_end] = feat.detach().float().cpu().numpy()
                if SAVE_N1:
                    n1_out[mb_beg:mb_end] = F.normalize(n_pred, dim=-1).detach().float().cpu().numpy()

            # 回填当前块到 memmap
            R = end - beg
            Fmm[beg:end, :, :] = feats_out.reshape(R, M, D_MODEL).astype(DTYPE_FEAT)
            Ymm[beg:end, :, :] = np.concatenate(Y_list, axis=0).reshape(R, M, 3).astype(DTYPE_Y)
            if SAVE_N1 and n1_out is not None:
                Nmm[beg:end, :, :] = n1_out.reshape(R, M, 3).astype(np.float16)

            pbar.update(R)

    Fmm.flush(); Ymm.flush()
    if SAVE_N1 and Nmm is not None:
        Nmm.flush()
    pbar.close()

    # 写 meta
    meta_out = {
        "mesh_name": os.path.basename(mesh_dir.rstrip(os.sep)),
        "d_model": D_MODEL,
        "lsd_r_size": int(meta['lsd_r_size']),
        "lsd_t_size": int(meta['lsd_t_size']),
        "N": 1 + int(meta['lsd_r_size']) * int(meta['lsd_t_size']),
        "K": int(K),
        "M": int(M),
        "dtype_feat": str(DTYPE_FEAT),
        "dtype_y": str(DTYPE_Y),
        "save_n1": SAVE_N1,
        "best_model_path": os.path.abspath(BEST_MODEL_PATH),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
    }
    with open(os.path.join(out_dir, "meta.json"), "w") as f:
        json.dump(meta_out, f, indent=2)

def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    meta = _load_meta(_default_meta_path(DATASET_ROOT))
    # 构建/加载 S1
    model = build_face_model(meta)

    mesh_dirs = _list_mesh_dirs(DATASET_ROOT)
    if not mesh_dirs:
        raise RuntimeError(f"No mesh dirs with lsd.npy & gt.npy under '{DATASET_ROOT}'")
    print(f"[meta] r={meta['lsd_r_size']}  t={meta['lsd_t_size']}  N={1+int(meta['lsd_r_size'])*int(meta['lsd_t_size'])}")
    print(f"[plan] meshes={len(mesh_dirs)}  out_root={OUT_ROOT}")
    print(f"[ckpt] {BEST_MODEL_PATH}")

    for mdir in mesh_dirs:
        mesh_name = os.path.basename(mdir.rstrip(os.sep))
        out_dir = os.path.join(OUT_ROOT, mesh_name)
        if os.path.exists(os.path.join(out_dir, "meta.json")):
            print(f"[skip] {mesh_name} already generated.")
            continue
        print(f"[gen ] {mesh_name}")
        process_one_mesh(mdir, PATCH_ROOT, out_dir, model, meta)
        print(f"[done] {mesh_name}")

if __name__ == "__main__":
    main()
