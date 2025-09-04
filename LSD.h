#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <unordered_map>
#include <unordered_set>
#include <stdio.h>
#include <map>
#include <queue>
#include <thread>
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <cstdint>
#include <cmath>

const int thread_number = 8;
extern std::thread td[thread_number];
// const int mt_flag = 0;
// const int lsdsize = 80;
const int lsd_r_size = 40;
const int lsd_t_size = 40;
const int sampling_size = lsd_r_size * lsd_t_size + 1;
const int ringnum = 4;
const int patch_num = 240;
// extern int supmat[lsdsize][lsdsize][3];

enum FaceNeighborType
{
	kVertexBased,
	kEdgeBased,
	kRadiusBased
};
enum DenoiseType
{
	kLocal,
	kGlobal
};

struct MyTraits : OpenMesh::DefaultTraits
{
	// Let Point and Normal be a vector of doubles
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	// The default 1D texture coordinate type is float.
	typedef double TexCoord1D;
	// The default 2D texture coordinate type is OpenMesh::Vec2f.
	typedef OpenMesh::Vec2d TexCoord2D;
	// The default 3D texture coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3d TexCoord3D;

	// enable standart properties
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
	HalfedgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::PrevHalfedge);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
	EdgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;

struct ring
{
	std::vector<int> facelist;
	std::vector<int> totalring[4];
};
struct line
{
	TriMesh::Point v1;
	TriMesh::Point v2;
};

struct SampleDirection
{
	float x;
	float y;
	float theta;
	float radius;
	float r2;
};

extern std::vector<SampleDirection> local_sample;

void getFaceNormal(TriMesh &mesh, std::vector<TriMesh::Normal> &normals);

void getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid);

double getSigmaS(double multiple, std::vector<TriMesh::Point> &centroid, TriMesh &mesh);

void makeRing(TriMesh &mesh, std::vector<ring> &ringlist, int ringnum);

bool CalculateLineLineIntersection(TriMesh::Point &line1Point1, TriMesh::Point &line1Point2,
								   TriMesh::Point &line2Point1, TriMesh::Point &line2Point2, TriMesh::Point &resultSegmentPoint, TriMesh::Normal &nownormal);

// void gsupmat();

TriMesh::Normal getAveNormal(
	const ring &cur_ring,
	const std::vector<TriMesh::Normal> &noisy_normals,
	int current_flag,
	const std::vector<int> &flagz);

TriMesh::Normal getPolarAxis(TriMesh &mesh, int face_index, const std::vector<TriMesh::Point> &face_centroid);

int samplingNormal(
	TriMesh &mesh,
	int index,
	const std::vector<TriMesh::Point> &face_centroid,
	const std::vector<TriMesh::Normal> &noisy_normals,
	std::vector<line> &halfedgeset,
	double sigma_s,
	std::vector<SampleDirection> &local_sample,
	float *outputmat);

std::vector<int> getPatch(TriMesh &mesh, int index, const int n_faces);

void markBoundaryFaces(TriMesh &mesh, std::vector<int> &flagz);

void generateLocalSamplingOrder(std::vector<SampleDirection> &local_sample);

// PoissonDiskCenterGen — direct (online) center generation by ring-distance Poisson disk
// - No IOU, no post-filtering. We forbid centers within r_block ring-distance.
// - Works during Gdata generation or denoising; O(#centers × (ring frontier)).
//
// Usage (pseudocode):
//   // You provide bfs_rings(center, max_rings, out) and get_patch(center, M, out)
//   PoissonDiskCenterGen gen(nfaces, /*K_patch=*/K, /*alpha=*/0.6, /*r_block_override=*/-1,
//                            /*shuffle=*/true, /*seed=*/42);
//   gen.set_bfs_rings([&](int c, int R, std::vector<int>& out){ bfs_rings(c, R, out); });
//   gen.set_get_patch([&](int c, int M, std::vector<int>& out){ get_patch_faces(c, M, out); });
//   auto centers = gen.run();
//   auto st = gen.stats(/*check_cover_with_patch=*/true, /*M=*/M);
//