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

const int thread_number = 8;
extern std::thread td[thread_number];
const int mt_flag = 1;  
const int lsdsize = 80;
extern float *outputcache;
extern int supmat[lsdsize][lsdsize][3];

enum FaceNeighborType { kVertexBased, kEdgeBased, kRadiusBased };
enum DenoiseType { kLocal, kGlobal };

struct MyTraits : OpenMesh::DefaultTraits
{
	// Let Point and Normal be a vector of doubles
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	// The default 1D texture coordinate type is float.
	typedef double  TexCoord1D;
	// The default 2D texture coordinate type is OpenMesh::Vec2f.
	typedef OpenMesh::Vec2d  TexCoord2D;
	// The default 3D texture coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3d  TexCoord3D;

	//enable standart properties
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


void getFaceNormal(TriMesh& mesh, std::vector<TriMesh::Normal>& normals);

void getFaceCentroid(TriMesh& mesh, std::vector<TriMesh::Point>& centroid);

double getSigmaS(double multiple, std::vector<TriMesh::Point>& centroid, TriMesh& mesh);

void makeRing(TriMesh &mesh, std::vector<ring> &ringlist, int ringnum);

bool CalculateLineLineIntersection(TriMesh::Point& line1Point1, TriMesh::Point& line1Point2,
	TriMesh::Point& line2Point1, TriMesh::Point& line2Point2, TriMesh::Point& resultSegmentPoint, TriMesh::Normal& nownormal);

void gsupmet();