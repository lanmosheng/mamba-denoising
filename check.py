#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, struct, ast, argparse

EXIT_OK           = 0
EXIT_FILE_MISSING = 5
EXIT_PARSE_ERROR  = 6
EXIT_SHAPE_BAD    = 2
EXIT_SIZE_SHORT   = 3
EXIT_SIZE_OVER    = 4

def read_npy_header(path):
    with open(path, 'rb') as f:
        magic = f.read(6)
        if magic != b'\x93NUMPY':
            raise ValueError("Not a .npy file (bad magic)")
        major, minor = struct.unpack('BB', f.read(2))
        if major == 1:
            hlen = struct.unpack('<H', f.read(2))[0]
        elif major in (2, 3):
            hlen = struct.unpack('<I', f.read(4))[0]
        else:
            raise ValueError(f"Unsupported .npy version: {major}.{minor}")
        header = f.read(hlen).decode('latin1')
        try:
            hdr = ast.literal_eval(header)
        except Exception as e:
            raise ValueError(f"Header parse failed: {e}")
        header_bytes = 6 + 2 + (2 if major == 1 else 4) + hlen
        return header_bytes, hdr

def human(n):
    # 简单的人类可读
    for unit in ['B','KB','MB','GB','TB']:
        if n < 1024.0:
            return f"{n:.2f}{unit}"
        n /= 1024.0
    return f"{n:.2f}PB"

def check_lsd(path, N_expected=1001, strict=True, verbose=True):
    if not os.path.isfile(path):
        if verbose:
            print(f"[error] file not found: {path}")
        return EXIT_FILE_MISSING

    try:
        header_bytes, hdr = read_npy_header(path)
    except Exception as e:
        print(f"[error] header read/parse: {e}")
        return EXIT_PARSE_ERROR

    shape = hdr.get('shape', None)
    descr = hdr.get('descr', None)
    if shape is None or descr is None:
        print(f"[error] header missing keys: have keys {list(hdr.keys())}")
        return EXIT_PARSE_ERROR

    try:
        itemsize = __import__('numpy').dtype(descr).itemsize
    except Exception:
        # 兜底：大多数情况 '<f4'
        itemsize = 4

    # 形状校验：(-1, 1001, 3)
    ok_shape = (isinstance(shape, tuple) and len(shape) == 3 and shape[1] == N_expected and shape[2] == 3)
    if not ok_shape:
        print(f"[shape] BAD: shape={shape} (expect (-1,{N_expected},3))")
        if strict:
            return EXIT_SHAPE_BAD
    else:
        print(f"[shape] OK : shape={shape}")

    nF_header = int(shape[0]) if (isinstance(shape, tuple) and len(shape) >= 1 and isinstance(shape[0], (int,)) and shape[0] >= 0) else None

    row_bytes = N_expected * 3 * itemsize
    filesize  = os.path.getsize(path)
    data_real = max(0, filesize - header_bytes)  # 实际数据区字节
    data_need = (nF_header * row_bytes) if (nF_header is not None) else None

    print(f"[file ] size={filesize} ({human(filesize)}), header_bytes={header_bytes}, row_bytes={row_bytes}")

    # 推算“按实际文件大小”已经完整写入的行数 / 下一行已写的 float 数
    full_rows   = data_real // row_bytes
    rem_bytes   = data_real %  row_bytes
    rem_floats  = rem_bytes // itemsize

    print(f"[rows ] by_size: full_rows={full_rows}, rem_floats_in_next_row={rem_floats}")

    if data_need is None:
        # 没有 nF（极少见），无法比较 size 是否足够
        print("[warn ] header has no concrete nF; cannot compare size vs need")
        return EXIT_OK

    diff = data_real - data_need
    if diff == 0:
        print("[size ] OK : data bytes match header declaration exactly")
        return EXIT_OK
    elif diff < 0:
        print(f"[size ] SHORT by {-diff} bytes ({-diff // itemsize} floats, {-diff // row_bytes} full rows)")
        if strict:
            return EXIT_SIZE_SHORT
        return EXIT_OK
    else:
        print(f"[size ] OVER by {diff} bytes ({diff // itemsize} floats, {diff // row_bytes} full rows)  <-- unexpected")
        return EXIT_SIZE_OVER

def main():
    ap = argparse.ArgumentParser(description="Check lsd.npy header/size consistency (after each block write).")
    ap.add_argument("npy_path", help="Path to lsd.npy")
    ap.add_argument("--N", type=int, default=1001, help="Expected second-dim (default: 1001)")
    ap.add_argument("--strict", action="store_true", help="Nonzero exit on shape/size mismatch")
    args = ap.parse_args()
    code = check_lsd(args.npy_path, N_expected=args.N, strict=args.strict, verbose=True)
    sys.exit(code)

if __name__ == "__main__":
    main()
