#include "LSD.h"
#include <fstream>
#include <cnpy.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#endif

static const size_t LSD_BLOCK_FACES = 1024; // 每次攒 1024 个面再写

// std::vector<float> outputcache;
std::vector<SampleDirection> local_sample;
struct pid
{
	int index;
	int meshindex;
	int count;
	pid()
	{
		index = 0;
		meshindex = 0;
		count = 0;
	}
	pid(int a, int b, int c)
	{
		index = a;
		meshindex = b;
		count = c;
	}
};
int genGt = 1;
std::vector<TriMesh> meshlist;

std::vector<pid> thread_p[thread_number];

std::vector<double> sigma_s_list;
std::vector<std::vector<ring>> ringlist_list;
std::vector<TriMesh> noisemeshlist;
std::vector<std::vector<line>> halfedgeset_list;
std::vector<std::vector<TriMesh::Normal>> noisy_normals_list;
std::vector<std::vector<TriMesh::Point>> face_centroid_list;
std::vector<std::vector<TriMesh::Normal>> filtered_normals_list;
std::vector<std::vector<int>> flagz_list;

int gLSD(int index, TriMesh &noisemesh, float outputmat[sampling_size * 3], float groundtruth[3],
		 double sigma_s,
		 std::vector<ring> &ringlist,
		 std::vector<TriMesh::Normal> &filtered_normals,
		 std::vector<line> &halfedgeset,
		 std::vector<TriMesh::Normal> &noisy_normals,
		 std::vector<TriMesh::Point> &face_centroid,
		 std::vector<int> &flagz)
{

	// obtain polar axis
	if (genGt)
	{
		Eigen::Vector3d gtnormal(filtered_normals[index].data()[0], filtered_normals[index].data()[1], filtered_normals[index].data()[2]);
		gtnormal.normalize();
		groundtruth[0] = (float)gtnormal[0];
		groundtruth[1] = (float)gtnormal[1];
		groundtruth[2] = (float)gtnormal[2];
	}

	// generate LSD
	int err = samplingNormal(noisemesh, index, face_centroid, noisy_normals, halfedgeset, sigma_s, local_sample, outputmat);
	return err;
}

int preprocessing(
	TriMesh &mesh,		// GT mesh
	TriMesh &noisemesh, // noisy mesh
	std::vector<ring> &ringlist,
	std::vector<TriMesh::Normal> &noisy_normals,
	std::vector<TriMesh::Normal> &filtered_normals,
	std::vector<TriMesh::Point> &face_centroid,
	std::vector<line> &halfedgeset,
	std::vector<int> &flagz,
	double &sigma_s,
	int nom)
{

	ringlist.resize(noisemesh.n_faces());
	noisy_normals.resize(noisemesh.n_faces());
	face_centroid.resize(noisemesh.n_faces());
	halfedgeset.resize(noisemesh.n_halfedges());

	for (TriMesh::HalfedgeIter it = noisemesh.halfedges_begin(); it != noisemesh.halfedges_end(); ++it)
	{
		halfedgeset[it->idx()].v1 = noisemesh.point(noisemesh.from_vertex_handle(*it));
		halfedgeset[it->idx()].v2 = noisemesh.point(noisemesh.to_vertex_handle(*it));
	}

	makeRing(noisemesh, ringlist, 3);
	getFaceNormal(noisemesh, noisy_normals);
	getFaceCentroid(noisemesh, face_centroid);
	sigma_s = getSigmaS(2, face_centroid, noisemesh);
	markBoundaryFaces(noisemesh, flagz);

	if (genGt == 1)
	{
		filtered_normals.resize(mesh.n_faces());
		getFaceNormal(mesh, filtered_normals);
	}

	return 0;
}

int mkfolder(std::string outputname)
{
	std::string dir = outputname; // 形如 "dataset/文件夹名/"
	while (!dir.empty() && (dir.back() == '/' || dir.back() == '\\'))
		dir.pop_back();

	int rc;
#ifdef _WIN32
	rc = _mkdir(dir.c_str());
#else
	rc = mkdir(dir.c_str(), 0755);
#endif
	if (rc != 0 && errno != EEXIST)
	{
		perror(("mkdir failed: " + dir).c_str());
		return 0;
	}
	return 1;
}

void generateFileBlock(const std::string &outdir,
					   const float *lsd_block, size_t K,
					   const float *gt_block, size_t Kgt,
					   bool first_block)
{
	const char *mode = first_block ? "w" : "a";
	std::string lsd_path = outdir + "/lsd.npy";
	// 一次写 K 行 (K, N, 3)
	cnpy::npy_save(lsd_path, lsd_block,
				   std::vector<size_t>{K, (size_t)sampling_size, (size_t)3},
				   mode);
	if (genGt == 1)
	{
		std::string gt_path = outdir + "/gt.npy";
		// 一次写 K 行 (K, 3)
		cnpy::npy_save(gt_path, gt_block,
					   std::vector<size_t>{Kgt, (size_t)3},
					   mode);
	}
}

// void generateFile(const std::string &outdir, float *lsdcache, float *gtcache, bool first)
// {
// 	const char *mode = first ? "w" : "a";
// 	std::string lsd_path = outdir + "/lsd.npy";
// 	cnpy::npy_save(lsd_path, lsdcache, std::vector<size_t>{1, sampling_size, 3}, mode); // (1,N,3)
// 	if (genGt == 1)
// 	{
// 		std::string gt_path = outdir + "/gt.npy";
// 		cnpy::npy_save(gt_path, gtcache, std::vector<size_t>{1, 3}, mode); // (1,3)
// 	}
// }
void generatePatchFile(const std::string &outdir, const std::vector<int> &patches, bool first)
{
	std::string patch_path = outdir + "/patch_faces.npy";
	const char *mode = first ? "w" : "a";
	std::vector<int32_t> row(patches.begin(), patches.end());
	cnpy::npy_save(patch_path, row.data(), std::vector<size_t>{1, (size_t)row.size()}, mode);
}
bool write_meta(const std::string &dir)
{
	std::string path = dir + "meta.json";
	std::ofstream ofs(path, std::ios::binary | std::ios::trunc);
	if (!ofs)
		return false;

	ofs << "{\n"
		<< "  \"lsd_r_size\": " << lsd_r_size << ",\n"
		<< "  \"lsd_t_size\": " << lsd_t_size << ",\n"
		<< "  \"ringnum\": " << ringnum << ",\n"
		<< "  \"patch_num\": " << patch_num << "\n"
		<< "}\n";
	ofs.close();
	return ofs.good();
}
void face_neighbors_by_edge(TriMesh &mesh, int f, std::vector<int> &out)
{
	out.clear();
	for (TriMesh::FaceHalfedgeIter fhe = mesh.fh_begin(TriMesh::FaceHandle(f));
		 fhe.is_valid(); ++fhe)
	{
		auto oppF = mesh.opposite_face_handle(*fhe);
		if (oppF.is_valid())
			out.push_back(oppF.idx());
	}
}
// B) 以“边相邻”为基础的环距 BFS：返回 [0..R] 内（闭包）的面集合
static void bfs_rings_edge(TriMesh &mesh, int center, int R, std::vector<int> &out)
{
	out.clear();
	const int nfaces = mesh.n_faces();
	std::vector<char> vis(nfaces, 0);
	std::vector<int> cur, nxt;
	cur.push_back(center);
	vis[center] = 1;
	out.push_back(center);

	for (int r = 0; r < R; ++r)
	{
		nxt.clear();
		for (int f : cur)
		{
			std::vector<int> nbs;
			face_neighbors_by_edge(mesh, f, nbs);
			for (int nb : nbs)
			{
				if (!vis[nb])
				{
					vis[nb] = 1;
					nxt.push_back(nb);
					out.push_back(nb);
				}
			}
		}
		if (nxt.empty())
			break;
		cur.swap(nxt);
	}
}

int main(int argc, char *argv[])
{
	int profile_num = 0;
	int numberofmesh = 0;

	srand(0);

	FILE *profile;
	if (argc == 2)
	{
		profile = fopen(argv[1], "r");
	}
	else
	{
		printf("profile error\n");
		return 0;
	}

	fscanf(profile, "%d", &genGt);

	fscanf(profile, "%d", &numberofmesh);
	meshlist.resize(numberofmesh);
	noisemeshlist.resize(numberofmesh);
	sigma_s_list.resize(numberofmesh);
	ringlist_list.resize(numberofmesh);
	filtered_normals_list.resize(numberofmesh);
	halfedgeset_list.resize(numberofmesh);
	noisy_normals_list.resize(numberofmesh);
	face_centroid_list.resize(numberofmesh);
	flagz_list.resize(numberofmesh);

	std::string outputfile;
	std::string patchfile;
	std::string mesh_n[numberofmesh * 2 + 1];

	// read noisy meshes
	for (int nom = 0; nom < numberofmesh; nom++)
	{
		char buff[30];
		fscanf(profile, "%s", buff);
		mesh_n[nom] = buff;
		if (!OpenMesh::IO::read_mesh(noisemeshlist[nom], mesh_n[nom]))
		{
			printf("read %s data error", mesh_n[nom].c_str());
			return 0;
		}
	}

	if (genGt == 1)
	{
		// read ground truth meshes
		printf("read mesh\n");
		for (int nom = 0; nom < numberofmesh; nom++)
		{
			char buff[30];
			fscanf(profile, "%s", buff);
			mesh_n[nom + numberofmesh] = buff;
			if (!OpenMesh::IO::read_mesh(meshlist[nom], mesh_n[nom + numberofmesh]))
			{
				printf("read %s mesh error", mesh_n[nom + numberofmesh].c_str());
				return 0;
			}
			if (noisemeshlist[nom].n_faces() != meshlist[nom].n_faces())
			{
				printf("read %s data error, number of faces differ", mesh_n[nom].c_str());
				return 0;
			}
		}
	}
	printf("read mesh over\n");

	int px[3]; // parameters for gdata
	// 0,1,2: the index range of output files groups, range(10, 20, 2) = 10, 12, 14, 16, 18

	char outputfilebuff[30];
	fscanf(profile, "%s", outputfilebuff); // name and path of dataset
	outputfile = outputfilebuff;

	char patchfilebuff[30];
	fscanf(profile, "%s", patchfilebuff);
	patchfile = patchfilebuff;

	for (int i = 0; i < 3; i++)
		fscanf(profile, "%d", &px[i]);

	int gen_lsd, gen_patch;

	fscanf(profile, "%d%d", &gen_lsd, &gen_patch);
	// printf("%d %d %d\n", px[0], px[1], px[2]);
	for (int nom = px[0]; nom < px[1]; nom += px[2])
	{
		preprocessing(
			meshlist[nom],
			noisemeshlist[nom],
			ringlist_list[nom],
			noisy_normals_list[nom],
			filtered_normals_list[nom],
			face_centroid_list[nom],
			halfedgeset_list[nom],
			flagz_list[nom],
			sigma_s_list[nom],
			nom);
	}

	printf("Write Meta\n");
	write_meta(outputfile);

	generateLocalSamplingOrder(local_sample);
	for (int k0 = px[0]; k0 < px[1]; k0 += px[2])
	{
		printf("Processing %s\n", mesh_n[k0].c_str());

		std::string name = mesh_n[k0];
		if (name.rfind("strain/", 0) == 0)
			name.erase(0, 7);
		auto pos = name.find_last_of('.'); // 找最后一个点
		if (pos != std::string::npos)
			name.erase(pos);

		int nfaces = noisemeshlist[k0].n_faces();
		if (gen_lsd)
		{
			std::string outputname = outputfile + name;
			mkfolder(outputname);
			printf("Generate LSD\n");

			float *lsd_block = new float[LSD_BLOCK_FACES * sampling_size * 3];
			float *gt_block = new float[LSD_BLOCK_FACES * 3];
			int cur_in_block = 0;
			bool first_block = true;
			for (int index = 0; index < nfaces; index++)
			{

				float *lsd_ptr = lsd_block + cur_in_block * (sampling_size * 3);
				float *gt_ptr = gt_block + cur_in_block * 3;

				if (gLSD(index, noisemeshlist[k0], lsd_ptr, gt_ptr, sigma_s_list[k0], ringlist_list[k0], filtered_normals_list[k0], halfedgeset_list[k0], noisy_normals_list[k0], face_centroid_list[k0], flagz_list[k0]) == -4)
				{
					printf("LSD Error %s faceindex %d\n", mesh_n[k0].c_str(), index);
					exit(1);
				}

				cur_in_block++;
				if (cur_in_block == LSD_BLOCK_FACES)
				{
					generateFileBlock(outputname, lsd_block, cur_in_block, gt_block, cur_in_block, first_block);
					first_block = false;
					cur_in_block = 0;
				}
			}

			if (cur_in_block > 0)
			{
				generateFileBlock(outputname,
								  lsd_block, cur_in_block,
								  gt_block, cur_in_block,
								  first_block);
			}

			delete[] lsd_block;
			delete[] gt_block;
		}
		if (gen_patch)
		{
			printf("Generate Patch (Poisson-Disk centers, edge-adj)\n");
			std::string patchname = patchfile + name;
			mkfolder(patchname);

			const int nfaces_local = nfaces;							 // 或 meshlist[k0].n_faces()
			const int K_patch = /* 你构建 patch 所用的“环层半径 K” */ 6; // 不确定就随便先给，或用 r_block_override
			const int M_expected = patch_num;

			// 1) 配置 Poisson 生成器（按边相邻）
			PoissonDiskCenterGen gen(
				nfaces_local,
				K_patch,
				0.5, // 阻隔系数：越大中心越稀疏；0.5~0.7 常用
				3,	 // 若不知道 K，可改成具体环距 r_block（如 5）
				true,
				42);

			// 2) 回调：环距 BFS 用“边相邻”，patch 仍用你现有 getPatch
			gen.set_bfs_rings([&](int c, int R, std::vector<int> &out)
							  {
								  bfs_rings_edge(noisemeshlist[k0], c, R, out); // ← 用上面的实现
							  });
			gen.set_get_patch([&](int c, int /*M*/, std::vector<int> &out)
							  {
								  out = getPatch(noisemeshlist[k0], c, nfaces_local); // 你已有的接口
							  });

			// 3) 在线选中心（不做 IOU、不过滤）
			const int max_centers = -1; // 可设上限避免爆时长；-1=不限
			const std::vector<int> &centers = gen.run(max_centers);

			// 4) （可选）补洞，确保 coverage≈1.0（薄片/孔洞时有用）
			// gen.fill_gaps_if_any(/*M=*/M_expected, /*r_block_fill=*/-1);

			// 5) 打印统计
			if (M_expected > 0)
			{
				auto st = gen.stats(/*check_cover_with_patch=*/true, /*M=*/M_expected);
				printf("[centers] r_block=%d | count=%zu | coverage=%.3f (%d/%d) | avg_cov=%.2f | max_cov=%d\n",
					   gen.r_block_used(), centers.size(), st.coverage_rate, st.covered_faces, nfaces_local,
					   st.avg_cov, st.max_cov);
			}
			else
			{
				printf("[centers] r_block=%d | count=%zu\n", gen.r_block_used(), centers.size());
			}

			// 6) 仅为“选中中心”写 patch
			for (size_t i = 0; i < centers.size(); ++i)
			{
				int c = centers[i];
				std::vector<int> patches = getPatch(noisemeshlist[k0], c, nfaces_local);
				generatePatchFile(patchname, patches, /*is_first=*/i == 0);
			}
		}
	}

	printf("Gdata Over!");
	return 0;
}
