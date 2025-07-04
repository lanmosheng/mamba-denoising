#include"LSD.h"

int supmat[lsdsize][lsdsize][3];
std::thread td[thread_number];
float *outputcache;
struct pid
{
	int index;
	int count;
	pid()
	{
		index = 0;
		count = 0;
	}
	pid(int a, int c)
	{
		index = a;
		count = c;
	}
};
std::vector<pid> thread_p[thread_number];
const int filesize = 10000;
double sigma_s = 0;
std::vector<ring> ringlist;
TriMesh noisemesh;
std::vector<line> halfedgeset;
std::vector<TriMesh::Normal> noisy_normals;
std::vector<TriMesh::Point> face_centroid;
std::vector<TriMesh::Normal> filtered_normals;
std::vector<Eigen::Matrix3d> msave;
std::vector<int> errorflag;
std::vector<int> flagz;
std::vector<FILE*> filepo;

int gLSD(int index, float outputmat[lsdsize*lsdsize * 3])
{
	TriMesh::Normal a1(0, 0, 0);

	//obtain n*
	for (int ii = 0; ii < ringlist[index].totalring[1].size(); ii++)
	{

		if (flagz[index] < 0)
		{
			if (flagz[ringlist[index].totalring[1][ii]] < 0)
				a1 += noisy_normals[ringlist[index].totalring[1][ii]];
		}
		else
		{
			a1 += noisy_normals[ringlist[index].totalring[1][ii]];
		}

	}
	a1.normalize();

	//obtain polar axis
	TriMesh::Point startpoint(0, 0, 0);
	int cc = 0;
	for (TriMesh::FaceVertexIter it = noisemesh.fv_begin(TriMesh::FaceHandle(index)); cc <= 1; cc++, it++)
	{
		startpoint += noisemesh.point(*it);
	}
	startpoint /= 2;
	TriMesh::Normal startnormal = startpoint - face_centroid[index];
	startnormal.normalize();


	//obtain rotation matrix and inverse rotation matrix
	Eigen::Matrix3d d2(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2]), Eigen::Vector3d(1, 0, 0)));

	Eigen::Matrix3d d2r(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2])));

	msave[index] = d2r;
	//generate LSD

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

			//generate the sampling points
			int halfedgenum = -1;
			int endflag = 0;
			OpenMesh::FaceHandle nowface(index);

			if (supmat[i][j][0] == 0 && supmat[i][j][1] == 0)
			{

				Eigen::Vector3d temp5(noisy_normals[nowface.idx()].data()[0], noisy_normals[nowface.idx()].data()[1],
					noisy_normals[nowface.idx()].data()[2]);
				temp5 = d2 * temp5;
				outputmat[i * lsdsize * 3 + j * 3] = (float)temp5[0];
				outputmat[i * lsdsize * 3 + j * 3 + 1] = (float)temp5[1];
				outputmat[i * lsdsize * 3 + j * 3 + 2] = (float)temp5[2];
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
					for (TriMesh::FaceHalfedgeIter it = noisemesh.fh_begin(nowface); it != noisemesh.fh_end(nowface); it++)
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
									outputmat[i * lsdsize * 3 + j * 3] = (float)temp5[0];
									outputmat[i * lsdsize * 3 + j * 3 + 1] = (float)temp5[1];
									outputmat[i * lsdsize * 3 + j * 3 + 2] = (float)temp5[2];

									endflag = -1;
									break;

								}

								//otherwise the geodesics should be extended to next face

								nowpoint = nextpoint;

								halfedgenum = noisemesh.opposite_halfedge_handle(*it).idx();
								OpenMesh::FaceHandle nextface = noisemesh.face_handle(noisemesh.opposite_halfedge_handle(*it));

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
	return 0;
}
void updateVertexPosition(TriMesh &mesh, std::vector<TriMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary)
{
	std::vector<TriMesh::Point> new_points(mesh.n_vertices());

	std::vector<TriMesh::Point> centroid;

	for (int iter = 0; iter < iteration_number; iter++)
	{
		getFaceCentroid(mesh, centroid);
		for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			TriMesh::Point p = mesh.point(*v_it);
			if (fixed_boundary && mesh.is_boundary(*v_it))
			{
				new_points.at(v_it->idx()) = p;
			}
			else
			{
				double face_num = 0.0;
				TriMesh::Point temp_point(0.0, 0.0, 0.0);
				for (TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++)
				{
					TriMesh::Normal temp_normal = filtered_normals[vf_it->idx()];
					TriMesh::Point temp_centroid = centroid[vf_it->idx()];
					temp_point += temp_normal * (temp_normal | (temp_centroid - p));
					face_num++;
				}
				p += temp_point / face_num;

				new_points.at(v_it->idx()) = p;
			}
		}

		for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
			mesh.set_point(*v_it, new_points[v_it->idx()]);
	}
}

int preprocessing()
{
	ringlist.resize(noisemesh.n_faces());
	noisy_normals.resize(noisemesh.n_faces());
	face_centroid.resize(noisemesh.n_faces());
	filtered_normals.resize(noisemesh.n_faces());
	errorflag.resize(noisemesh.n_faces());
	msave.resize(noisemesh.n_faces());
	halfedgeset.resize(noisemesh.n_halfedges());
	for (TriMesh::HalfedgeIter it = noisemesh.halfedges_begin(); it != noisemesh.halfedges_end(); it++)
	{
		halfedgeset[(*it).idx()].v1 = noisemesh.point(noisemesh.from_vertex_handle(*it));
		halfedgeset[(*it).idx()].v2 = noisemesh.point(noisemesh.to_vertex_handle(*it));
	}
	makeRing(noisemesh, ringlist, 3);
	getFaceNormal(noisemesh, noisy_normals);
	getFaceCentroid(noisemesh, face_centroid);
	sigma_s = getSigmaS(2, face_centroid, noisemesh) / 8;
	flagz.resize(noisemesh.n_faces());
	for (TriMesh::FaceIter v_it = noisemesh.faces_begin(); v_it != noisemesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		flagz[index] = 0;
	}
	for (TriMesh::FaceIter v_it = noisemesh.faces_begin(); v_it != noisemesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		int count = 0;
		for (TriMesh::FaceFaceIter ff_it = noisemesh.ff_begin(TriMesh::FaceHandle(*v_it)); ff_it.is_valid(); ff_it++)
			count++;
		if (count <= 2)
		{
			flagz[index] = -1;
		}
	} //find the faces on the boundary and mark them with -1
	for (TriMesh::FaceIter v_it = noisemesh.faces_begin(); v_it != noisemesh.faces_end(); v_it++)
	{

		int index = v_it->idx();
		if (flagz[index] <= -1)
			continue;
		for (TriMesh::FaceVertexIter fv_it = noisemesh.fv_begin(TriMesh::FaceHandle(index)); fv_it.is_valid(); fv_it++)
		{

			for (TriMesh::VertexFaceIter vf_it = noisemesh.vf_begin(*fv_it); vf_it.is_valid(); vf_it++)
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


	for (TriMesh::FaceIter v_it = noisemesh.faces_begin(); v_it != noisemesh.faces_end(); v_it++)
	{
		int index = v_it->idx();
		if (flagz[index] == -2)
		{
			flagz[index] = -1;
		}
	}
	return 0;
}

int goutputfile(int nface, int ncase)
{
	int nfile = nface%ncase ? nface / ncase + 1 : nface / ncase;

	filepo.resize(nfile);
	FILE *namelist = fopen("list.txt", "w");
	fprintf(namelist, "%d\n", nfile);
	for (int i = 0; i < nfile; i++)
	{
		std::string  temps = std::to_string(i) + ".bin";
		filepo[i] = fopen(temps.c_str(), "wb");
		fprintf(namelist, "%s\n", temps.c_str());
	}
	fclose(namelist);
	return 0;
}
int closeall(std::vector<FILE*> &filepo)
{
	for (int i = 0; i < filepo.size(); i++)
	{
		fclose(filepo[i]);
	}
	return 0;
}

// generate the name of the denoised mesh
std::string gofn(std::string input, int ifv)
{
	int len = input.length();

	std::string temp = input.substr(0, len - 4);
	if (len < 7)
	{
		std::string temp3 = "_00.off";
		int itn = ifv + 1;
		temp3[1] = itn / 10 + '0';
		temp3[2] = itn % 10 + '0';
		std::string output = temp + temp3;
		return output;
	}
	std::string temp2 = input.substr(len - 4 - 3, 3);

	if (temp2[0] == '_'&&temp2[1] >= '0'&&temp2[1] <= '9'&&temp2[2] >= '0'&&temp2[2] <= '9')
	{

		int itn = (temp2[1] - '0') * 10 + (temp2[2] - '0');
		itn += ifv + 1;
		temp2[1] = itn / 10 + '0';
		temp2[2] = itn % 10 + '0';

		std::string output = input.substr(0, len - 7) + temp2 + ".off";
		return output;
	}
	else
	{

		std::string temp3 = "_00.off";
		int itn = ifv + 1;
		temp3[1] = itn / 10 + '0';
		temp3[2] = itn % 10 + '0';
		std::string output = temp + temp3;
		return output;
	}
}
void threadprocess(int p)
{
	for (int i = 0; i < thread_p[p].size(); i++)
	{
		int index = thread_p[p][i].index;
		int count = thread_p[p][i].count;

		if (gLSD(index, outputcache + count*lsdsize*lsdsize * 3) == -4)
			errorflag[index] = 1;
		else
			errorflag[index] = 0;

	}

}
int main(int argc, char* argv[])
{
	int profile_num = 0;
	int numberofmesh = 0;
	int numberofmodel = 0;
	int ifn, ivn;
	char modelpath[5][100];
	srand(0);

	FILE* profile;
	if (argc == 2)
	{
		profile = fopen(argv[1], "r");
	}
	else
	{
		printf("profile error\n");
		return 0;
	}



	fscanf(profile, "%d", &numberofmodel);
	for (int i = 0; i < numberofmodel; i++)
	{
		fscanf(profile, "%s", modelpath[i]);
	}

	fscanf(profile, "%d %d", &ifn, &ivn);
	if (ifn > numberofmodel || ivn <= 0)
	{
		printf("iter number error");
		return 0;
	}

	fscanf(profile, "%d", &numberofmesh);

	gsupmet();

	//read ground truth meshes 
	printf("read mesh\n");
	noisemesh.clean();
	outputcache = new float[filesize * lsdsize*lsdsize * 3];
	for (int nom = 0; nom < numberofmesh; nom++)
	{

		int totalfilenumber = 0;
		char mesh_n[100];
		fscanf(profile, "%s", mesh_n);
		printf("processing: ");
		printf("%s\n", mesh_n);
		if (!OpenMesh::IO::read_mesh(noisemesh, mesh_n))
		{
			printf("data error");
			return 0;
		}

		for (int iter = 0; iter < ifn; iter++)
		{
			double sigma_s = 0;
			ringlist.clear();
			halfedgeset.clear();
			noisy_normals.clear();
			face_centroid.clear();
			filtered_normals.clear();
			msave.clear();
			errorflag.clear();
			flagz.clear();
			filepo.clear();
			preprocessing();
			goutputfile(noisemesh.n_faces(), filesize);
			int count = 0;
			int fcount = 0;

			memset(outputcache, 0, filesize * lsdsize*lsdsize * 3 * sizeof(float));
			for (int k1 = 0; k1 < thread_number; k1++)
				thread_p[k1].clear();
			//write the LSD of all the faces to files
			int n_faces = noisemesh.n_faces();
			for (int iterf = 0; iterf < n_faces;)
			{
				thread_p[count%thread_number].push_back(pid(iterf, count));
				count++;
				iterf++;
				if (count == filesize || iterf == n_faces)
				{
					for (int k2 = 0; k2 < thread_number; k2++)
					{
						td[k2] = std::thread(threadprocess, k2);
					}
					for (int k2 = 0; k2 < thread_number; k2++)
					{
						td[k2].join();
					}
					fwrite(outputcache, sizeof(float), count * lsdsize*lsdsize * 3, filepo[fcount]);

					count = 0;
					fcount++;

					memset(outputcache, 0, filesize * lsdsize*lsdsize * 3 * sizeof(float));
					for (int k1 = 0; k1 < thread_number; k1++)
						thread_p[k1].clear();
				}


			}

			closeall(filepo);

			//call python to compute the normalized denoised normals
			std::string pycmd = "python denoising.py ";
			pycmd += std::string(modelpath[iter]) + " list.txt";

			system(pycmd.c_str());

			// read the denoised normals and denormalize them
			float *nomralcache = new float[noisemesh.n_faces() * 3];
			FILE *nf = fopen("normal.bin", "rb");
			fread(nomralcache, sizeof(float), noisemesh.n_faces() * 3, nf);
			fclose(nf);

			for (int iterf = 0; iterf<noisemesh.n_faces(); iterf++)
			{
				if (errorflag[iterf] == 0)
				{
					Eigen::Vector3d tt(nomralcache[iterf * 3], nomralcache[iterf * 3 + 1], nomralcache[iterf * 3 + 2]);
					tt = msave[iterf] * tt;

					filtered_normals[iterf] = TriMesh::Point(tt.data()[0], tt.data()[1], tt.data()[2]);
					filtered_normals[iterf].normalized();
				}

				else
				{
					filtered_normals[iterf] = noisy_normals[iterf];
				}
			}
			delete nomralcache;

			
			updateVertexPosition(noisemesh, filtered_normals, ivn, false);
			std::string outfilename = gofn(mesh_n, iter);
			OpenMesh::IO::write_mesh(noisemesh, outfilename);
		}
		noisemesh.clean();
	}
	delete outputcache;
	return 0;
}

