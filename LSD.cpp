#include"LSD.h"

void getFaceNormal(TriMesh& mesh, std::vector<TriMesh::Normal>& normals)
{
	mesh.request_face_normals();
	mesh.update_face_normals();

	normals.resize(mesh.n_faces());
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Normal n = mesh.normal(*f_it);
		normals[f_it->idx()] = n;
	}
}

void getFaceCentroid(TriMesh& mesh, std::vector<TriMesh::Point>& centroid)
{
	centroid.resize(mesh.n_faces(), TriMesh::Point(0.0, 0.0, 0.0));
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point pt = mesh.calc_face_centroid(*f_it);
		centroid[(*f_it).idx()] = pt;
	}
}

double getSigmaS(double multiple, std::vector<TriMesh::Point>& centroid, TriMesh& mesh)
{
	double sigma_s = 0.0, num = 0.0;
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point fi = centroid[f_it->idx()];
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
		{
			TriMesh::Point fj = centroid[ff_it->idx()];
			sigma_s += (fj - fi).length();
			num++;
		}
	}
	return sigma_s * multiple / num;
}

void makeRing(TriMesh &mesh, std::vector<ring> &ringlist, int ringnum)
{
	std::set<int> neighbor_face_index;

	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		int t = f_it->idx();
		ringlist[t].facelist.clear();

		for (int i = 0; i <ringnum; i++)
		{
			neighbor_face_index.clear();
			ringlist[t].facelist.push_back(f_it->idx());
			for (std::vector<int>::iterator iter = ringlist[t].facelist.begin(); iter != ringlist[t].facelist.end(); ++iter)
			{
				for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(TriMesh::FaceHandle(*iter)); fv_it.is_valid(); fv_it++)
				{
					for (TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*fv_it); vf_it.is_valid(); vf_it++)
						neighbor_face_index.insert(vf_it->idx());
				}

			}
			ringlist[t].facelist.clear();
			for (std::set<int>::iterator iter = neighbor_face_index.begin(); iter != neighbor_face_index.end(); ++iter)
			{
				ringlist[t].facelist.push_back(*iter);
			}

			for (int j = 0; j < ringlist[t].facelist.size(); j++)
				ringlist[t].totalring[i].push_back(ringlist[t].facelist[j]);



		}
	}
}

bool CalculateLineLineIntersection(TriMesh::Point& line1Point1, TriMesh::Point& line1Point2,
	TriMesh::Point& line2Point1, TriMesh::Point& line2Point2, TriMesh::Point& resultSegmentPoint, TriMesh::Normal& nownormal)
{
	TriMesh::Point p1 = line1Point1;
	TriMesh::Point p2 = line1Point2;
	TriMesh::Point p3 = line2Point1;
	TriMesh::Point p4 = line2Point2;
	TriMesh::Point p13 = p1 - p3;
	TriMesh::Point p43 = p4 - p3;

	if (p43.length() < 1e-8) {
		return false;
	}
	TriMesh::Point p21 = p2 - p1;
	if (p21.length() < 1e-8) {
		return false;
	}

	double d1343 = p13.data()[0] * (double)p43.data()[0] + (double)p13.data()[1] * p43.data()[1] + (double)p13.data()[2] * p43.data()[2];
	double d4321 = p43.data()[0] * (double)p21.data()[0] + (double)p43.data()[1] * p21.data()[1] + (double)p43.data()[2] * p21.data()[2];
	double d1321 = p13.data()[0] * (double)p21.data()[0] + (double)p13.data()[1] * p21.data()[1] + (double)p13.data()[2] * p21.data()[2];
	double d4343 = p43.data()[0] * (double)p43.data()[0] + (double)p43.data()[1] * p43.data()[1] + (double)p43.data()[2] * p43.data()[2];
	double d2121 = p21.data()[0] * (double)p21.data()[0] + (double)p21.data()[1] * p21.data()[1] + (double)p21.data()[2] * p21.data()[2];

	double denom = d2121 * d4343 - d4321 * d4321;
	if (denom == 0)
	{
		return false;
	}
	double numer = d1343 * d4321 - d1321 * d4343;

	double mua = numer / denom;

	double mub = (d1343 + d4321 * (mua)) / d4343;
	TriMesh::Point resultSegmentPoint1(
		(p1.data()[0] + mua * p21.data()[0]),
		(p1.data()[1] + mua * p21.data()[1]),
		(p1.data()[2] + mua * p21.data()[2]));
	TriMesh::Point resultSegmentPoint2(
		(p3.data()[0] + mub * p43.data()[0]),
		(p3.data()[1] + mub * p43.data()[1]),
		(p3.data()[2] + mub * p43.data()[2]));
	if ((resultSegmentPoint2 - resultSegmentPoint1).length() < 1e-6 && mua >= -1e-6 && mua <= 1 + 1e-6 && mub >= 0)
	{

		if (mua > 1 - 1e-6)
		{
			resultSegmentPoint1 = TriMesh::Point(
				(p1.data()[0] + 0.9999 * p21.data()[0]),
				(p1.data()[1] + 0.9999 * p21.data()[1]),
				(p1.data()[2] + 0.9999 * p21.data()[2]));
			nownormal = (resultSegmentPoint1 - line2Point1);
			nownormal.normalize();
		}
		if (mua < 1e-6)
		{
			resultSegmentPoint1 = TriMesh::Point(
				(p1.data()[0] + 0.0001 * p21.data()[0]),
				(p1.data()[1] + 0.0001 * p21.data()[1]),
				(p1.data()[2] + 0.0001 * p21.data()[2]));
			nownormal = (resultSegmentPoint1 - line2Point1);
			nownormal.normalize();
		}
		resultSegmentPoint = resultSegmentPoint1;
		return true;
	}
	else
		return false;
}

void gsupmet()
{
	memset(supmat, 0, sizeof(supmat));
	//generate i, j and i^2+j^2 
	for (int i = 0; i < lsdsize; i++)
		for (int j = 0; j < lsdsize; j++)
		{

		int zi = i - lsdsize / 2;

		int zj = j - lsdsize / 2;
		supmat[i][j][2] = zi * zi + zj * zj;
		supmat[i][j][0] = zi;
		supmat[i][j][1] = zj;

		}
	return;
}