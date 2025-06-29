import os
import open3d as o3d
import numpy as np

def load_mesh(file_path):
    mesh = o3d.io.read_triangle_mesh(file_path)
    mesh.compute_vertex_normals()
    return mesh

def chamfer_distance(pc1, pc2):
    # pc1, pc2 are open3d.geometry.PointCloud
    pcd_tree = o3d.geometry.KDTreeFlann(pc2)
    dist1 = []
    for point in pc1.points:
        [_, idx, d] = pcd_tree.search_knn_vector_3d(point, 1)
        dist1.append(np.sqrt(d[0]))

    pcd_tree = o3d.geometry.KDTreeFlann(pc1)
    dist2 = []
    for point in pc2.points:
        [_, idx, d] = pcd_tree.search_knn_vector_3d(point, 1)
        dist2.append(np.sqrt(d[0]))

    return np.mean(dist1), np.mean(dist2)

def main(folder):
    results = []

    for file in os.listdir(folder):
        if file.endswith("_01.off"):
            denoised_path = os.path.join(folder, file)
            # 推测对应的 gt 是类似 elephant.obj（去掉 _n1_01）
            basename = file.replace("_n1_01.off", "").replace("_01.off", "")
            gt_path = os.path.join(folder, basename + ".obj")

            if not os.path.exists(gt_path):
                print(f"Ground truth not found for {file}, expected: {gt_path}")
                continue

            mesh_denoised = load_mesh(denoised_path)
            mesh_gt = load_mesh(gt_path)

            # 采样点云
            pcd_denoised = mesh_denoised.sample_points_uniformly(number_of_points=10000)
            pcd_gt = mesh_gt.sample_points_uniformly(number_of_points=10000)

            cd1, cd2 = chamfer_distance(pcd_denoised, pcd_gt)
            avg_cd = (cd1 + cd2) / 2
            print(f"{file} → Chamfer Distances: GT→D {cd1:.6f}, D→GT {cd2:.6f}, Avg: {avg_cd:.6f}")
            results.append((file, cd1, cd2, avg_cd))

    print("\n=== Summary ===")
    for f, d1, d2, avg in results:
        print(f"{f}: Avg CD = {avg:.6f}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python evaluate_denoising.py <path_to_folder>")
    else:
        main(sys.argv[1])
