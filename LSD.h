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
const int lsd_r_size = 25;
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

struct PDStats
{
	double coverage_rate = 0.0; // fraction of faces covered by >=1 selected patch (needs get_patch to compute)
	int covered_faces = 0;
	double avg_cov = 0.0; // mean coverage per face
	int max_cov = 0;	  // max coverage per face
};

class PoissonDiskCenterGen
{
public:
	using BfsRingsFn = std::function<void(int /*center*/, int /*max_rings*/, std::vector<int> & /*out*/)>;
	using GetPatchFn = std::function<void(int /*center*/, int /*M*/, std::vector<int> & /*out*/)>;

	PoissonDiskCenterGen(int nfaces, int K_patch, double alpha = 0.6,
						 int r_block_override = -1, bool shuffle = true, uint64_t seed = 42)
		: nfaces_(nfaces), K_(K_patch), alpha_(alpha), r_block_override_(r_block_override),
		  shuffle_(shuffle), rng_(seed)
	{
		order_.resize(nfaces_);
		for (int i = 0; i < nfaces_; ++i)
			order_[i] = i;
		if (shuffle_)
			std::shuffle(order_.begin(), order_.end(), rng_);
		blocked_.assign(nfaces_, 0u);
	}

	void reseed(uint64_t seed)
	{
		rng_.seed(seed);
		if (shuffle_)
			std::shuffle(order_.begin(), order_.end(), rng_);
	}

	void set_bfs_rings(BfsRingsFn fn) { bfs_rings_ = std::move(fn); }
	void set_get_patch(GetPatchFn fn) { get_patch_ = std::move(fn); }

	// primary generation
	const std::vector<int> &run(int max_centers = -1)
	{
		ensure_bfs();
		const int r_block = (r_block_override_ > 0 ? r_block_override_ : std::max(1, (int)std::ceil(alpha_ * K_)));
		centers_.clear();

		std::vector<int> ringbuf;
		ringbuf.reserve(8 * K_); // heuristic

		for (int c : order_)
		{
			if (blocked_[c])
				continue; // already in forbidden zone
			// Accept c as a new center
			centers_.push_back(c);
			// Block a disk of radius r_block in ring metric
			ringbuf.clear();
			bfs_rings_(c, r_block, ringbuf);
			for (int v : ringbuf)
				if ((unsigned)v < (unsigned)nfaces_)
					blocked_[v] = 1u;
			if (max_centers > 0 && (int)centers_.size() >= max_centers)
				break;
		}
		return centers_;
	}

	// Optional gap-filling pass: ensure every face is covered by >=1 patch (using get_patch)
	// Strategy: greedy add a center at any still-uncovered face (respecting a smaller block, r_block_fill)
	// NOTE: call after run(). This keeps Poisson quality while plugging holes if your mesh has thin parts.
	void fill_gaps_if_any(int M, int r_block_fill = -1)
	{
		ensure_get_patch();
		ensure_bfs();
		if (r_block_fill <= 0)
			r_block_fill = std::max(1, (int)std::ceil(0.4 * K_));

		// Build coverage from existing centers
		std::vector<int> coverage(nfaces_, 0);
		std::vector<int> buf;
		buf.reserve(M);
		for (int c : centers_)
		{
			buf.clear();
			get_patch_(c, M, buf);
			for (int f : buf)
				if ((unsigned)f < (unsigned)nfaces_)
					++coverage[f];
		}

		// Scan faces; if uncovered, force-select as new center and block a smaller disk
		std::vector<int> ringbuf;
		ringbuf.reserve(8 * K_);
		for (int f : order_)
		{
			if (coverage[f] > 0)
				continue;
			centers_.push_back(f);
			buf.clear();
			get_patch_(f, M, buf);
			for (int v : buf)
				if ((unsigned)v < (unsigned)nfaces_)
					++coverage[v];
			ringbuf.clear();
			bfs_rings_(f, r_block_fill, ringbuf);
			for (int v : ringbuf)
				if ((unsigned)v < (unsigned)nfaces_)
					blocked_[v] = 1u;
		}
	}

	PDStats stats(bool check_cover_with_patch, int M) const
	{
		PDStats s{};
		if (!check_cover_with_patch)
			return s;
		if (!get_patch_)
			return s;
		std::vector<int> cov(nfaces_, 0), buf;
		buf.reserve(M);
		for (int c : centers_)
		{
			buf.clear();
			get_patch_(c, M, buf);
			for (int f : buf)
				if ((unsigned)f < (unsigned)nfaces_)
					++cov[f];
		}
		long long sum = 0;
		int covered = 0, maxcov = 0;
		for (int v = 0; v < nfaces_; ++v)
		{
			int x = cov[v];
			sum += x;
			if (x > 0)
				++covered;
			if (x > maxcov)
				maxcov = x;
		}
		s.covered_faces = covered;
		s.coverage_rate = (nfaces_ > 0 ? double(covered) / double(nfaces_) : 0.0);
		s.avg_cov = (nfaces_ > 0 ? double(sum) / double(nfaces_) : 0.0);
		s.max_cov = maxcov;
		return s;
	}

	const std::vector<int> &centers() const { return centers_; }
	int r_block_used() const { return (r_block_override_ > 0 ? r_block_override_ : std::max(1, (int)std::ceil(alpha_ * K_))); }

private:
	void ensure_bfs() const
	{
		if (!bfs_rings_)
			throw std::runtime_error("PoissonDiskCenterGen: bfs_rings not set");
	}
	void ensure_get_patch() const
	{
		if (!get_patch_)
			throw std::runtime_error("PoissonDiskCenterGen: get_patch not set");
	}

	int nfaces_ = 0;			// #faces in mesh
	int K_ = 1;					// BFS ring radius used to form a patch (approximate)
	double alpha_ = 0.6;		// block radius = ceil(alpha * K)
	int r_block_override_ = -1; // if >0, use this instead of alpha*K
	bool shuffle_ = true;
	mutable std::mt19937_64 rng_;

	std::vector<uint8_t> blocked_;
	std::vector<int> order_;
	std::vector<int> centers_;

	BfsRingsFn bfs_rings_ = nullptr;
	GetPatchFn get_patch_ = nullptr; // optional, only for stats / gap-fill
};
