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

void gsupmat()
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


TriMesh::Normal getAveNormal(
    const ring& cur_ring, 
    const std::vector<TriMesh::Normal> &noisy_normals, 
    int current_flag,
    const std::vector<int>& flagz) 
{
    TriMesh::Normal a1(0,0,0);
    for(int ii = 0; ii < cur_ring.totalring[1].size(); ii++){
        int neighbor_face = cur_ring.totalring[1][ii];
        if(current_flag < 0){
            if(flagz[neighbor_face] < 0)
                a1 += noisy_normals[cur_ring.totalring[1][ii]];
        }
        else{
            a1 += noisy_normals[neighbor_face];
        }
    }
    a1.normalize();
    return a1;
}


TriMesh::Normal getPolarAxis(TriMesh& mesh, int face_index, const std::vector<TriMesh::Point>& face_centroid)
{
	TriMesh::Point startpoint(0, 0, 0);
    int cc = 0;

    for (TriMesh::FaceVertexIter it = mesh.fv_begin(TriMesh::FaceHandle(face_index)); cc <= 1; cc++, it++) {
        startpoint += mesh.point(*it);
    }

    startpoint /= 2.0;
    TriMesh::Normal startnormal = startpoint - face_centroid[face_index];
    startnormal.normalize();

    return startnormal;
}

int samplingNormal(
	TriMesh& mesh,
    int index,
	const Eigen::Matrix3d& d2,
    const TriMesh::Normal& startnormal,
    const std::vector<TriMesh::Point>& face_centroid,
    const std::vector<TriMesh::Normal>& noisy_normals,
    std::vector<line>& halfedgeset,
    double sigma_s,
    int lsdsize,
    float* outputmat)
{
	for (int i = 0; i < lsdsize; i++)
	{
		for (int j = 0; j < lsdsize; j++)
		{

			double glength = sigma_s * sqrt(supmat[i][j][2]); //total length of geodesics
			double clength = 0; // current length of geodesics
			int centreface = index;

			TriMesh::Point nowpoint = face_centroid[index];
			TriMesh::Normal nownormal = startnormal;//direction of geodesics

			//compute the direction of geodesics
			if (supmat[i][j][0] == 0)
			{
				if (supmat[i][j][1] > 0)
				{
					Eigen::AngleAxisd rotation_vector(0, Eigen::Vector3d(noisy_normals[index].data()[0],
						noisy_normals[index].data()[1],
						noisy_normals[index].data()[2]));

					Eigen::Vector3d temp3(nownormal.data()[0], nownormal.data()[1], nownormal.data()[2]);
					temp3 = rotation_vector.matrix() * temp3;
					nownormal.data()[0] = temp3[0];
					nownormal.data()[1] = temp3[1];
					nownormal.data()[2] = temp3[2];
				}
				else
				{
					Eigen::AngleAxisd rotation_vector(M_PI, Eigen::Vector3d(noisy_normals[index].data()[0],
						noisy_normals[index].data()[1],
						noisy_normals[index].data()[2]));

					Eigen::Vector3d temp3(nownormal.data()[0], nownormal.data()[1], nownormal.data()[2]);
					temp3 = rotation_vector.matrix() * temp3;
					nownormal.data()[0] = temp3[0];
					nownormal.data()[1] = temp3[1];
					nownormal.data()[2] = temp3[2];
				}

			}
			else
			{
				Eigen::AngleAxisd rotation_vector(atan2(supmat[i][j][0], supmat[i][j][1]), Eigen::Vector3d(noisy_normals[index].data()[0],
					noisy_normals[index].data()[1],
					noisy_normals[index].data()[2]));
				Eigen::Vector3d temp3(nownormal.data()[0], nownormal.data()[1], nownormal.data()[2]);
				temp3 = rotation_vector.matrix() * temp3;
				nownormal.data()[0] = temp3[0];
				nownormal.data()[1] = temp3[1];
				nownormal.data()[2] = temp3[2];
			}
			nownormal.normalize();

			//generate the normal of sampling points
			int halfedgenum = -1;
			int endflag = 0;
			OpenMesh::FaceHandle nowface(index);

			if (supmat[i][j][0] == 0 && supmat[i][j][1] == 0)
			{

				Eigen::Vector3d temp5(noisy_normals[nowface.idx()].data()[0], noisy_normals[nowface.idx()].data()[1],
					noisy_normals[nowface.idx()].data()[2]);
				temp5 = d2 * temp5;
				outputmat[i * lsdsize*3 + j * 3] = (float)temp5[0];
				outputmat[i * lsdsize*3 + j * 3+ 1] = (float)temp5[1];
				outputmat[i * lsdsize*3 + j * 3+ 2] = (float)temp5[2];
			}
			else
			{
				std::unordered_set<int> visitedface;

				while (endflag == 0)
				{
					visitedface.insert(nowface.idx());
					int edgecount = 0;
					int goflag = 0;
					//find the edge that intersects the geodesic
					for (TriMesh::FaceHalfedgeIter it = mesh.fh_begin(nowface); it != mesh.fh_end(nowface); it++)
					{
						TriMesh::Point nextpoint, temppoint;
						int nowhalfedge = (*it).idx();


						if (nowhalfedge != halfedgenum)
						{
							edgecount++;
							auto temppoint = nowpoint + nownormal * sigma_s * 100;

							if (CalculateLineLineIntersection(halfedgeset[nowhalfedge].v1,
								halfedgeset[nowhalfedge].v2, nowpoint, temppoint, nextpoint, nownormal) == true)
							{

								//found the intersection
								goflag = 1;
								clength += (nextpoint - nowpoint).length();

								// if the current length of geodesics is longer than the total length
								// the sampling point is located at this face
								// save the face normal to outputmat
								if (clength >= glength)
								{
									Eigen::Vector3d temp5(noisy_normals[nowface.idx()].data()[0], noisy_normals[nowface.idx()].data()[1],
										noisy_normals[nowface.idx()].data()[2]);
									temp5 = d2 * temp5;
									temp5.normalize();
									outputmat[i * lsdsize*3 + j * 3] = (float)temp5[0];
									outputmat[i * lsdsize * 3 + j * 3 + 1] = (float)temp5[1];
									outputmat[i * lsdsize * 3 + j * 3 + 2] = (float)temp5[2];


									endflag = -1;
									break;

								}

								//otherwise the geodesics should be extended to next face

								nowpoint = nextpoint;

								halfedgenum = mesh.opposite_halfedge_handle(*it).idx();
								OpenMesh::FaceHandle nextface = mesh.face_handle(mesh.opposite_halfedge_handle(*it));

								if (nextface.idx() == -1)  //reach the boundray, stop 
								{
									endflag = -2;
									break;
								}
								if (visitedface.find(nextface.idx()) != visitedface.end()) // reach a visited face, stop
								{
									endflag = -3;
									break;
								}
								// rotate the geodesics to next face
								Eigen::Matrix3d d(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(noisy_normals[nowface.idx()].data()[0],
									noisy_normals[nowface.idx()].data()[1],
									noisy_normals[nowface.idx()].data()[2]), Eigen::Vector3d(noisy_normals[nextface.idx()].data()[0],
									noisy_normals[nextface.idx()].data()[1],
									noisy_normals[nextface.idx()].data()[2])));

								Eigen::Vector3d temp2(nownormal.data()[0], nownormal.data()[1], nownormal.data()[2]);
								temp2 = d * temp2;
								nownormal.data()[0] = temp2[0];
								nownormal.data()[1] = temp2[1];
								nownormal.data()[2] = temp2[2];
								nownormal.normalize();

								nowface = nextface;
								break;
							}

						}

						// can not find the edge, error
						if ((edgecount == 2 && goflag == 0 && halfedgenum != -1) || (edgecount == 3 && goflag == 0 && halfedgenum == -1))
						{
							endflag = -4;
							printf("error: %d %d %d|", index, i, j);
							return endflag;
						}
					}

				}
			}

		}
	}
}

std::vector<int> globalSampling(TriMesh& mesh, const std::vector<int>& flagz, const int n_faces){
	std::unordered_set<int> visited;
	std::vector<int> ret;
	
	auto bfs = [&](int start){
		std::queue<int> q;
		q.push(start);
		visited.insert(start);

		while(!q.empty()){
			int cur = q.front();
			q.pop();
			ret.push_back(cur);

			for(TriMesh::FaceFaceIter ff_it = mesh.ff_begin(TriMesh::FaceHandle(cur)); ff_it.is_valid(); ff_it++){
				int nxt = ff_it -> idx();
				if(!visited.count(nxt)){
					visited.insert(nxt);
					q.push(nxt);
				}
			}
		}
	};


	for(int i = 0; i < n_faces; i++){
		if(visited.count(i)) continue;
		if(flagz[i] < 0) continue;
		bfs(i);
	}

	for(int i = 0; i < n_faces; i++){
		if(visited.count(i)) continue;
		bfs(i);
	}

	if (ret.size() != n_faces) {
		std::cerr << "Warning: not all faces included in sampling!" << std::endl;
	}

	return ret;
}

void markBoundaryFaces(TriMesh& mesh, std::vector<int>& flagz) {
    int n_faces = mesh.n_faces();
    flagz.resize(n_faces);
	for (TriMesh::FaceIter v_it = mesh.faces_begin(); v_it != mesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		flagz[index] = 0;
	}
	for (TriMesh::FaceIter v_it = mesh.faces_begin(); v_it != mesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		int count = 0;
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_begin(TriMesh::FaceHandle(*v_it)); ff_it.is_valid(); ff_it++)
			count++;
		if (count <= 2)
		{
			flagz[index] = -1;
		}
	} //find the faces on the boundary and mark them with -1
	for (TriMesh::FaceIter v_it = mesh.faces_begin(); v_it != mesh.faces_end(); v_it++)
	{

		int index = v_it->idx();
		if (flagz[index] <= -1)
			continue;
		for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(TriMesh::FaceHandle(index)); fv_it.is_valid(); fv_it++)
		{

			for (TriMesh::VertexFaceIter vf_it = mesh.vf_begin(*fv_it); vf_it.is_valid(); vf_it++)
			{
				if (flagz[vf_it->idx()] == -1)
				{
					flagz[index] = -2;
					break;
				}
			}
			if (flagz[index] == -2)
				break;
		}

	} //find the faces that their neighbours are on the boundary, and mark them with -2


	for (TriMesh::FaceIter v_it = mesh.faces_begin(); v_it != mesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		if (flagz[index] == -2)
		{
			flagz[index] = -1;
		}
	}
}
