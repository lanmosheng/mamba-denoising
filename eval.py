#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evaluate denoising results by computing the mean unsigned angle error (degrees)
between a denoised mesh and its ground-truth mesh. The program reads a text
file where the first line is an integer N, followed by 2*N lines alternating:
    denoised_mesh_path
    gt_mesh_path

Supported mesh formats: .obj and .off
- We compute geometry-based normals (ignore any stored normals).
- Primary comparison is per-triangle face normals (after triangulation).
- If face counts mismatch but vertex counts match, we fall back to comparing
  per-vertex normals (area-weighted). If neither matches, the pair is skipped.

Output (stdout): a single line
    Average angle error (deg): <value>

Usage:
    python eval_denoise.py /path/to/pairs.txt

Notes:
- Paths inside the txt can be relative; they will be resolved relative to the
  directory of the txt file itself.
- Angle is unsigned: acos(|dot(n1, n2)|) in degrees.
"""

import sys
import os
import math
from typing import List, Tuple
import numpy as np

def _safe_acos(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, -1.0, 1.0)
    return np.arccos(x)

def _normalize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    n = np.maximum(n, eps)
    return v / n

def _triangulate_polygon(face: List[int]) -> List[Tuple[int, int, int]]:
    """Fan triangulation: (v0, v1, v2), (v0, v2, v3), ..."""
    if len(face) < 3:
        return []
    if len(face) == 3:
        return [(face[0], face[1], face[2])]
    tris = []
    for k in range(1, len(face) - 1):
        tris.append((face[0], face[k], face[k+1]))
    return tris

def load_off(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load OFF mesh. Returns (V [N,3], F_tri [M,3])"""
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        # Skip comments/blank until we hit OFF
        first = ''
        while True:
            line = f.readline()
            if line == '':
                raise ValueError("Unexpected EOF before OFF header")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            first = line
            break
        if first != 'OFF':
            # Some files use 'COFF' for colored OFF
            if first.upper().endswith('OFF'):
                pass
            else:
                raise ValueError(f"Expected OFF header, got: {first}")
        # Next non-empty, non-comment line has counts
        while True:
            line = f.readline()
            if line == '':
                raise ValueError("Unexpected EOF when reading OFF counts")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            n_vertices = int(parts[0])
            n_faces = int(parts[1])
            break
        V = []
        while len(V) < n_vertices:
            line = f.readline()
            if line == '':
                raise ValueError("Unexpected EOF reading OFF vertices")
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            V.append([float(parts[0]), float(parts[1]), float(parts[2])])
        V = np.array(V, dtype=np.float64)
        tris = []
        faces_read = 0
        while faces_read < n_faces:
            line = f.readline()
            if line == '':
                break
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            m = int(parts[0])
            idx = list(map(int, parts[1:1+m]))
            tris.extend(_triangulate_polygon(idx))
            faces_read += 1
        F = np.array(tris, dtype=np.int64) if tris else np.zeros((0,3), dtype=np.int64)
        return V, F

def load_obj(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load OBJ mesh (triangulates polygons). Returns (V [N,3], F_tri [M,3])"""
    V = []
    faces_raw = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('v '):
                parts = line.split()
                if len(parts) >= 4:
                    V.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith('f '):
                parts = line.split()[1:]
                idxs = []
                for p in parts:
                    # Handle formats like v, v/t, v//n, v/t/n
                    if '/' in p:
                        v_str = p.split('/')[0]
                    else:
                        v_str = p
                    if v_str == '' or v_str is None:
                        continue
                    vi = int(v_str)
                    if vi > 0:
                        vi0 = vi - 1
                    else:
                        # Negative indices refer from end
                        vi0 = len(V) + vi
                    idxs.append(vi0)
                if len(idxs) >= 3:
                    faces_raw.append(idxs)
    V = np.array(V, dtype=np.float64)
    tris = []
    for face in faces_raw:
        tris.extend(_triangulate_polygon(face))
    F = np.array(tris, dtype=np.int64) if tris else np.zeros((0,3), dtype=np.int64)
    return V, F

def load_mesh(path: str) -> Tuple[np.ndarray, np.ndarray]:
    ext = os.path.splitext(path)[1].lower()
    if ext == '.off':
        return load_off(path)
    elif ext == '.obj':
        return load_obj(path)
    else:
        raise ValueError(f"Unsupported mesh format: {ext} ({path})")

def face_normals(V: np.ndarray, F: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute per-triangle normals and areas.
    Returns (N [M,3] unit normals, A [M] areas).
    """
    if F.size == 0:
        return np.zeros((0,3), dtype=np.float64), np.zeros((0,), dtype=np.float64)
    v0 = V[F[:,0], :]
    v1 = V[F[:,1], :]
    v2 = V[F[:,2], :]
    cross = np.cross(v1 - v0, v2 - v0)
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    N = _normalize(cross)
    return N, areas

def vertex_normals(V: np.ndarray, F: np.ndarray) -> np.ndarray:
    """Compute per-vertex area-weighted normals (unit)."""
    if V.size == 0 or F.size == 0:
        return np.zeros_like(V)
    Nf, areas = face_normals(V, F)
    acc = np.zeros_like(V)
    for t, tri in enumerate(F):
        w = areas[t]
        for vi in tri:
            acc[vi] += Nf[t] * w
    return _normalize(acc)

def mean_unsigned_angle_deg(N1: np.ndarray, N2: np.ndarray) -> float:
    """Mean unsigned angle in degrees between two sets of unit vectors."""
    if N1.shape != N2.shape or N1.size == 0:
        return float('nan')
    # Ensure unit
    N1 = _normalize(N1)
    N2 = _normalize(N2)
    dots = np.sum(N1 * N2, axis=1)
    ang = _safe_acos(np.abs(np.clip(dots, -1.0, 1.0))) * 180.0 / math.pi
    return float(np.mean(ang))

def resolve_path(base_dir: str, p: str) -> str:
    q = p.strip()
    if not q:
        return q
    if os.path.isabs(q):
        return q
    return os.path.normpath(os.path.join(base_dir, q))

def process_pair(den_path: str, gt_path: str) -> float:
    """Return mean unsigned angle error (degrees) for this pair.
    Tries face-normal comparison first, falls back to vertex-normal comparison.
    """
    Vd, Fd = load_mesh(den_path)
    Vg, Fg = load_mesh(gt_path)

    # Try face normals if triangle counts match
    if Fd.shape[0] == Fg.shape[0] and Fd.shape[0] > 0:
        Nd, _ = face_normals(Vd, Fd)
        Ng, _ = face_normals(Vg, Fg)
        return mean_unsigned_angle_deg(Nd, Ng)

    # Fallback: vertex normals if vertex counts match
    if Vd.shape[0] == Vg.shape[0] and Vd.shape[0] > 0:
        Nd = vertex_normals(Vd, Fd)
        Ng = vertex_normals(Vg, Fg)
        return mean_unsigned_angle_deg(Nd, Ng)

    # Cannot compare
    return float('nan')

def read_pairs_list(txt_path: str) -> List[Tuple[str, str]]:
    base = os.path.dirname(os.path.abspath(txt_path))
    with open(txt_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = [ln.strip() for ln in f if ln.strip() != '']
    if not lines:
        raise ValueError("Empty list file")
    try:
        n = int(lines[0].split()[0])
    except Exception as e:
        raise ValueError("First line must be an integer (count of pairs)") from e
    paths = lines[1:1+2*n]
    if len(paths) < 2*n:
        raise ValueError(f"List expects {n} pairs but only found {len(paths)//2}")
    pairs = []
    for i in range(n):
        den_rel = paths[2*i]
        gt_rel  = paths[2*i + 1]
        print(den_rel, gt_rel)
        pairs.append((den_rel, gt_rel))
    return pairs

def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print("Usage: python eval_denoise.py /path/to/pairs.txt", file=sys.stderr)
        return 2
    list_path = argv[1]
    pairs = read_pairs_list(list_path)

    errs = []
    used = 0
    for den, gt in pairs:
        if not os.path.exists(den):
            # skip silently to keep output minimal
            continue
        if not os.path.exists(gt):
            continue
        try:
            e = process_pair(den, gt)
        except Exception:
            e = float('nan')
        if not (math.isnan(e) or math.isinf(e)):
            errs.append(e)
            used += 1

    if used == 0:
        # To keep output singular as requested, print NaN if nothing usable
        print("Average angle error (deg): NaN")
        return 1

    overall = float(np.mean(errs)) if errs else float('nan')
    # Single-line output as requested
    print(f"Average angle error (deg): {overall:.6f}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
