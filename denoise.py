#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, math, torch, json
import numpy as np
import torch.nn.functional as F

# ====== 路径与开关 ======
MESH_ROOT   = 'stest'
LSD_ROOT    = 'valset/test1001'
GT_ROOT     = 'valset/test1001'     # 与训练一致，用于评估
PATCH_ROOT  = 'valpatches'
OUTPUT_PATH = 'f_result'

# 优先选择本次训练产物；若不存在，会自动降级
MODEL_PATH_CANDIDATES = [
    'out/face_agg/stage3_best.pt',   # 联合微调最优
    'out/face_agg/stage2_best.pt',   # 仅Patch-Encoder最优
    'out/face_agg/stage1_best.pt',   # 仅Face-Encoder最优
]

# 顶点位置更新相关参数
POSITION_ITERS = 20     # 迭代次数（3~20可调）
FIXED_BOUNDARY = True   # 是否固定边界顶点

from collections import Counter

def _build_vertex_adjacency(Fidx, nV):
    """顶点 -> 相邻面 列表"""
    adj = [[] for _ in range(nV)]
    for f, (i0, i1, i2) in enumerate(Fidx):
        adj[i0].append(f); adj[i1].append(f); adj[i2].append(f)
    return adj

def _boundary_vertices(Fidx, nV):
    """通过“只被一个三角形使用”的边来判定边界顶点"""
    edges = []
    for i0, i1, i2 in Fidx:
        edges.append(tuple(sorted((i0, i1))))
        edges.append(tuple(sorted((i1, i2))))
        edges.append(tuple(sorted((i2, i0))))
    cnt = Counter(edges)
    bmask = np.zeros(nV, dtype=bool)
    for (a, b), c in cnt.items():
        if c == 1:  # 边只出现一次 => 边界
            bmask[a] = True; bmask[b] = True
    return bmask

def _row_norm(x, eps=1e-12):
    n = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(n, eps, None)

def eval_angle_metrics_deg(pred_face, gt_face, cover_mask=None, unsigned=True):
    """
    pred_face, gt_face: (nF,3)
    unsigned=True 时用 arccos(|dot|)
    返回: dict(count, mean, median, p90, p95, max)
    """
    pf = _row_norm(pred_face.astype(np.float64))
    gf = _row_norm(gt_face.astype(np.float64))

    if unsigned:
        dots = np.abs(np.sum(pf * gf, axis=1))
    else:
        dots = np.sum(pf * gf, axis=1)
    dots = np.clip(dots, -1.0, 1.0)
    errs = np.degrees(np.arccos(dots))  # (nF,)

    if cover_mask is not None:
        m = cover_mask.astype(bool)
        errs = errs[m]
    if errs.size == 0:
        return dict(count=0, mean=np.nan, median=np.nan, p90=np.nan, p95=np.nan, max=np.nan)

    return dict(
        count=int(errs.size),
        mean=float(np.mean(errs)),
        median=float(np.median(errs)),
        p90=float(np.percentile(errs, 90)),
        p95=float(np.percentile(errs, 95)),
        max=float(np.max(errs)),
    )

def update_vertex_positions(V, Fidx, face_normals, iters=5, fixed_boundary=True):
    """
    等价于你的 C++: updateVertexPosition
    V: (nV,3) 顶点坐标
    Fidx: (nF,3) 面顶点索引
    face_normals: (nF,3) 去噪后的面法线（单位向量）
    返回: (nV,3) 更新后的顶点坐标
    """
    V = V.astype(np.float32).copy()
    nV = V.shape[0]
    # 邻接与边界准备
    adj   = _build_vertex_adjacency(Fidx, nV)
    bmask = _boundary_vertices(Fidx, nV) if fixed_boundary else np.zeros(nV, dtype=bool)

    FN = face_normals.astype(np.float64)
    for _ in range(iters):
        centroids = (V[Fidx[:, 0]] + V[Fidx[:, 1]] + V[Fidx[:, 2]]) / 3.0
        centroids = centroids.astype(np.float64)

        delta = np.zeros((nV, 3), dtype=np.float64)
        count = np.zeros((nV,),   dtype=np.int32)

        for f, (i0, i1, i2) in enumerate(Fidx):
            n = FN[f]
            n_norm = np.linalg.norm(n)
            if n_norm > 0:
                n = n / n_norm
            c = centroids[f]
            for i in (i0, i1, i2):
                p = V[i].astype(np.float64)
                d = float(np.dot(n, (c - p)))   # 标量
                delta[i] += n * d
                count[i] += 1

        up_idx = np.where(~bmask & (count > 0))[0]
        V[up_idx] = (V[up_idx].astype(np.float64) + (delta[up_idx] / count[up_idx, None])).astype(np.float32)

    return V

# === 读取简易 OBJ（v / f） ===
def load_obj_simple(path):
    vs, fs = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            if line.startswith('v '):
                _, x, y, z = line.strip().split()[:4]
                vs.append([float(x), float(y), float(z)])
            elif line.startswith('f '):
                items = line.strip().split()[1:]
                idx = []
                for it in items:
                    v_str = it.split('/')[0]
                    if v_str == '':
                        continue
                    idx.append(int(v_str) - 1)  # 1-based -> 0-based
                if len(idx) == 3:
                    fs.append(idx)
                else:
                    for k in range(1, len(idx) - 1):
                        fs.append([idx[0], idx[k], idx[k + 1]])
    V = np.asarray(vs, dtype=np.float32)   # (nV,3)
    Fidx = np.asarray(fs, dtype=np.int32)  # (nF,3)
    return V, Fidx

# === 角度加权，把面法线聚合成顶点法线（导出 OBJ 用） ===
def vertex_normals_from_face_normals(V, F, NF):
    nV = V.shape[0]
    VN = np.zeros((nV, 3), dtype=np.float32)
    for (i0, i1, i2), n in zip(F, NF):
        v0, v1, v2 = V[i0], V[i1], V[i2]
        def angle(a, b, c):
            u = a - b; v = c - b
            u = u / (np.linalg.norm(u) + 1e-12)
            v = v / (np.linalg.norm(v) + 1e-12)
            cos = np.clip((u * v).sum(), -1.0, 1.0)
            return math.acos(cos)
        a0 = angle(v1, v0, v2)
        a1 = angle(v2, v1, v0)
        a2 = angle(v0, v2, v1)
        VN[i0] += a0 * n
        VN[i1] += a1 * n
        VN[i2] += a2 * n
    nn = np.linalg.norm(VN, axis=1, keepdims=True)
    VN = VN / np.clip(nn, 1e-12, None)
    return VN

# === 写 OBJ：vn 与 v 同索引，f 行用 v//vn ===
def write_obj_with_normals(path, V, F, VN):
    with open(path, "w") as f:
        for x, y, z in V:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for nx, ny, nz in VN:
            f.write(f"vn {nx:.6f} {ny:.6f} {nz:.6f}\n")
        for i0, i1, i2 in F:
            a = i0 + 1; b = i1 + 1; c = i2 + 1
            f.write(f"f {a}//{a} {b}//{b} {c}//{c}\n")

# === 几何面法线（用于与预测对齐朝向） ===
def face_normals_from_geometry(V, F):
    v0 = V[F[:, 0]]; v1 = V[F[:, 1]]; v2 = V[F[:, 2]]
    n  = np.cross(v1 - v0, v2 - v0)          # (nF,3)
    nn = np.linalg.norm(n, axis=1, keepdims=True)
    n  = n / np.clip(nn, 1e-12, None)
    return n.astype(np.float32)

# === Rodrigues 批量 R(a->ex)，与训练一致：锚旋到 +X ===
def rodrigues_to_ex(anchor, eps=1e-6, device="cpu", dtype=torch.float32):
    B = anchor.shape[0]
    a = F.normalize(anchor, dim=-1, eps=eps)
    b = torch.tensor([1.0, 0.0, 0.0], device=device, dtype=dtype).expand(B, 3)
    v = torch.cross(a, b, dim=-1)
    c = (a * b).sum(dim=-1)
    s = torch.linalg.norm(v, dim=-1)
    I = torch.eye(3, device=device, dtype=dtype).unsqueeze(0).expand(B, 3, 3)
    near0 = (s < 1e-12)

    alt_x = torch.tensor([1, 0, 0], device=device, dtype=dtype).expand_as(a)
    alt_y = torch.tensor([0, 1, 0], device=device, dtype=dtype).expand_as(a)
    use_alt_y = (a[:, 0].abs() > 0.9).unsqueeze(-1)
    alt = torch.where(use_alt_y, alt_y, alt_x)
    axis_anti = torch.cross(a, alt, dim=-1)
    axis_anti = F.normalize(axis_anti, dim=-1, eps=eps)
    zeros = torch.zeros(B, device=device, dtype=dtype)
    K180 = torch.stack([
        zeros,            -axis_anti[:, 2],  axis_anti[:, 1],
        axis_anti[:, 2],   zeros,           -axis_anti[:, 0],
       -axis_anti[:, 1],   axis_anti[:, 0],  zeros
    ], dim=1).reshape(B, 3, 3)
    R_anti = I + 2.0 * (K180 @ K180)

    axis = v / s.clamp_min(eps).unsqueeze(-1)
    K = torch.stack([
        zeros,       -axis[:, 2],  axis[:, 1],
        axis[:, 2],   zeros,      -axis[:, 0],
       -axis[:, 1],   axis[:, 0],  zeros
    ], dim=1).reshape(B, 3, 3)
    R_general = I + K * s.view(-1, 1, 1) + (K @ K) * (1.0 - c).view(-1, 1, 1)

    R = torch.where(near0.view(-1, 1, 1) & (c > 0).view(-1, 1, 1), I, R_general)
    R = torch.where(near0.view(-1, 1, 1) & (~(c > 0)).view(-1, 1, 1), R_anti, R)
    return R

# === 面聚合：把 (B,M,3) 预测 + face_idx(B,M) 聚合为 (nF,3) ===
def aggregate_by_face(pred_faces_all, face_idx_all, n_faces, return_counts=True):
    accum = np.zeros((n_faces, 3), dtype=np.float64)
    count = np.zeros((n_faces, 1), dtype=np.float64)
    for P, I in zip(pred_faces_all, face_idx_all):
        P = np.asarray(P, dtype=np.float64)
        I = np.asarray(I, dtype=np.int64)
        np.add.at(accum, I, P)
        np.add.at(count, I, 1.0)
    used = (count[:, 0] > 0)
    out = np.zeros((n_faces, 3), dtype=np.float32)
    out[used] = (accum[used] / np.clip(count[used], 1.0, None)).astype(np.float32)
    out /= np.clip(np.linalg.norm(out, axis=1, keepdims=True), 1e-12, None)
    if return_counts:
        return out, count[:, 0]
    return out

# === 选择可用的 MODEL_PATH（按候选顺序） ===
def _select_model_path():
    for p in MODEL_PATH_CANDIDATES:
        if os.path.isfile(p):
            return p
    # 兜底提示：老路径可能是旧架构（例如 ResNet 的），不推荐
    # 但为保持兼容性，如果你确实要用旧的，可在此返回旧路径
    raise FileNotFoundError(
        f"No model file found in candidates: {MODEL_PATH_CANDIDATES}. "
        f"Please set a correct checkpoint path."
    )

# === 模型载入：使用与你训练一致的构建方式，并修正去前缀逻辑 ===
def load_model(ckpt_path, device, get_model_fn):
    model = get_model_fn().to(device)
    sd = torch.load(ckpt_path, map_location=device)
    if isinstance(sd, dict) and 'model' in sd:
        sd = sd['model']
    # 正确处理 'module.' 前缀
    new_sd = {}
    for k, v in sd.items():
        if k.startswith('module.'):
            new_sd[k[7:]] = v
        else:
            new_sd[k] = v
    missing, unexpected = model.load_state_dict(new_sd, strict=False)
    print(f"[load] keys={len(new_sd)} missing={len(missing)} unexpected={len(unexpected)}")
    if missing:
        print("  missing (head):", list(missing)[:10])
    if unexpected:
        print("  unexpected(head):", list(unexpected)[:10])
    model.eval()
    return model

# === 用你项目里的 config.get_model 保持与训练完全一致 ===
def get_model_from_project():
    meta_path = os.path.join(LSD_ROOT, 'meta.json')
    with open(meta_path, 'r', encoding='utf-8') as f:
        meta = json.load(f)
    from train_utils import config as cfg
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = cfg.get_model(
        device,
        patch_num=meta['patch_num'],
        lsd_r_size=meta['lsd_r_size'],
        lsd_t_size=meta['lsd_t_size'],
        # 若你在 train.py 里覆写了这些默认值，请在此显式保持一致：
        # d_model=64, face_depth=4, patch_depth=4,
        # add_pos_emb=False, add_patch_pos=False,
        # d_state=16, d_conv=4, expand=2, residual_scale=0.5, p_drop=0.05,
    )
    return model

def mean_angle_deg(A, B, eps=1e-8):
    A = _row_norm(A.astype(np.float64))
    B = _row_norm(B.astype(np.float64))
    dot = np.sum(A * B, axis=1)
    dot = np.clip(dot, -1.0, 1.0)
    return float(np.degrees(np.arccos(dot)).mean())

# === 主逻辑 ===
def run():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model_path = _select_model_path()
    print(f"[info] Using model: {model_path}")
    model = load_model(model_path, device, get_model_from_project)

    mesh_names = sorted([d for d in os.listdir(LSD_ROOT) if os.path.isdir(os.path.join(LSD_ROOT, d))])
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    for name in mesh_names:
        obj_path   = os.path.join(MESH_ROOT,  f"{name}.obj")
        lsd_path   = os.path.join(LSD_ROOT,   name, "lsd.npy")
        patch_path = os.path.join(PATCH_ROOT, name, "patch_faces.npy")

        if not (os.path.isfile(obj_path) and os.path.isfile(lsd_path) and os.path.isfile(patch_path)):
            print(f"[skip] {name}: missing file(s)")
            continue

        V, Face = load_obj_simple(obj_path)
        nF = Face.shape[0]

        lsd = np.load(lsd_path).astype(np.float32)            # (nF, N, 3)
        patch_faces = np.load(patch_path).astype(np.int64)    # (P, M)
        P, M = patch_faces.shape
        N = lsd.shape[1]

        pred_chunks = []
        faceid_chunks = []

        with torch.no_grad():
            BATCH = 8
            for s in range(0, P, BATCH):
                e = min(P, s + BATCH)
                sel = patch_faces[s:e]                # (B, M)
                Xb = lsd[sel, :, :]                   # (B, M, N, 3)
                Xb = torch.from_numpy(Xb).to(device)
                # 前处理：与训练一致（旋到 +X 并归一化）
                Xb = F.normalize(Xb, dim=-1, eps=1e-6)

                anchor = Xb[:, 0, 0, :]               # (B, 3) patch中心面的 LSD[0]
                R = rodrigues_to_ex(anchor, device=device)
                Rt = R.transpose(1, 2)

                XT = Xb.view(Xb.shape[0], -1, 3)      # (B, M*N, 3)
                X_rot = torch.einsum('bij,bkj->bki', Rt, XT).view_as(Xb)

                out = model(X_rot)                    # TwoStageMamba => (n_hat, n1)
                y_rot = out[0] if isinstance(out, tuple) else out  # (B, M, 3)

                # 反旋回到全局
                yT = torch.einsum('bij,bkj->bki', R, y_rot)
                y_glb = F.normalize(yT, dim=-1, eps=1e-6).cpu().numpy()

                # 收集
                pred_chunks.extend([y_glb[i] for i in range(y_glb.shape[0])])
                faceid_chunks.extend([sel[i] for i in range(sel.shape[0])])

        # 面级聚合（与训练评估口径一致）
        pred_face, cover = aggregate_by_face(pred_chunks, faceid_chunks, nF)
        cover_mask = (cover > 0)

        # 评估（若 GT 可用）
        gt_path = os.path.join(GT_ROOT, name, "gt.npy")
        if os.path.isfile(gt_path):
            gt = np.load(gt_path).astype(np.float32)  # (nF,3)
            m = eval_angle_metrics_deg(pred_face, gt, cover_mask=cover_mask, unsigned=True)
            print(f"[{name}] GT unsigned deg: mean={m['mean']:.3f}  median={m['median']:.3f}  "
                  f"p90={m['p90']:.3f}  p95={m['p95']:.3f}  max={m['max']:.3f}  (used={m['count']}/{nF})")
        else:
            print(f"[{name}] GT not found, skip eval: {gt_path}")

        # 使预测法线与几何法线朝向一致（便于顶点更新/导出）
        n_orig = face_normals_from_geometry(V, Face)
        flip = (pred_face * n_orig).sum(axis=1) < 0.0
        pred_face[flip] = -pred_face[flip]

        # 顶点更新
        V_updated = update_vertex_positions(
            V, Face, pred_face,
            iters=POSITION_ITERS,
            fixed_boundary=FIXED_BOUNDARY
        )

        # 写出 OBJ（顶点法线为角度加权聚合）
        VN = vertex_normals_from_face_normals(V_updated, Face, pred_face)
        out_obj = os.path.join(OUTPUT_PATH, f"{name}_denoised.obj")
        write_obj_with_normals(out_obj, V_updated, Face, VN)

        print(f"[ok] {name} -> {out_obj}")

if __name__ == "__main__":
    run()
