#include"LSD.h"

int supmat[lsdsize][lsdsize][3];
std::thread td[thread_number];
float *outputcache;
float *gtcache;

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
std::vector<pid> thread_p[thread_number];


std::vector<TriMesh> meshlist;
std::vector<TriMesh> noisemeshlist;
std::vector<double> sigma_s;
std::vector<std::vector<ring>> ringlist;
std::vector<std::vector<TriMesh::Normal>> filtered_normals;
std::vector<std::vector<line>> halfedgeset;
std::vector<std::vector<TriMesh::Normal>> noisy_normals;
std::vector<std::vector<TriMesh::Point>> face_centroid;
std::vector<std::vector<int>> flagz;

int gLSD(int index, TriMesh &mesh2, float outputmat[lsdsize*lsdsize*3], float groundtruth[3],
	double sigma_s,
	std::vector<ring> &ringlist,
	std::vector<TriMesh::Normal> &filtered_normals,
	std::vector<line> &halfedgeset,
	std::vector<TriMesh::Normal> &noisy_normals,
	std::vector<TriMesh::Point> &face_centroid,
	std::vector<int> &flagz)
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
	for (TriMesh::FaceVertexIter it = mesh2.fv_begin(TriMesh::FaceHandle(index)); cc <= 1; cc++, it++)
	{
		startpoint += mesh2.point(*it);
	}
	startpoint /= 2;
	TriMesh::Normal startnormal = startpoint - face_centroid[index];
	startnormal.normalize();


	//obtain rotation matrix and rotated ground truth
	Eigen::Matrix3d d2(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2]), Eigen::Vector3d(1, 0, 0)));

	Eigen::Vector3d gtnormal(filtered_normals[index].data()[0], filtered_normals[index].data()[1], filtered_normals[index].data()[2]);
	gtnormal = d2 * gtnormal;
	gtnormal.normalize();

	groundtruth[0] = (float)gtnormal[0];
	groundtruth[1] = (float)gtnormal[1];
	groundtruth[2] = (float)gtnormal[2];

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
					for (TriMesh::FaceHalfedgeIter it = mesh2.fh_begin(nowface); it != mesh2.fh_end(nowface); it++)
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

								halfedgenum = mesh2.opposite_halfedge_handle(*it).idx();
								OpenMesh::FaceHandle nextface = mesh2.face_handle(mesh2.opposite_halfedge_handle(*it));

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

void threadprocess(int p)
{
	for (int i = 0; i < thread_p[p].size(); i++)
	{
		int index = thread_p[p][i].index;
		int meshidx = thread_p[p][i].meshindex;
		int count = thread_p[p][i].count;

		gLSD(index, noisemeshlist[meshidx], outputcache + count*lsdsize*lsdsize * 3, gtcache + count * 3, sigma_s[meshidx], ringlist[meshidx], filtered_normals[meshidx], halfedgeset[meshidx], noisy_normals[meshidx], face_centroid[meshidx], flagz[meshidx]);

	}

}
int main(int argc, char* argv[])
{
	int profile_num = 0;
	int numberofmesh = 0;
	
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
	fscanf(profile, "%d", &numberofmesh);
	meshlist.resize(numberofmesh);
	noisemeshlist.resize(numberofmesh);
	sigma_s.resize(numberofmesh);
	ringlist.resize(numberofmesh);
	filtered_normals.resize(numberofmesh);
	halfedgeset.resize(numberofmesh);
	noisy_normals.resize(numberofmesh);
	face_centroid.resize(numberofmesh);
	flagz.resize(numberofmesh);

	std::vector<std::pair<int, int>> traindata;

	char outputname[200];
	char outputflagname[200];
	char mesh_n[200];


	//read ground truth meshes 
	printf("read mesh\n");
	for (int nom = 0; nom < numberofmesh; nom++)
	{

		fscanf(profile, "%s", mesh_n);
		if (!OpenMesh::IO::read_mesh(meshlist[nom], mesh_n))
		{
			printf("data error");
			return 0;
		}

		ringlist[nom].resize(meshlist[nom].n_faces());
		noisy_normals[nom].resize(meshlist[nom].n_faces());
		face_centroid[nom].resize(meshlist[nom].n_faces());
		filtered_normals[nom].resize(meshlist[nom].n_faces());
	}
	//read noisy meshes 
	for (int nom = 0; nom < numberofmesh; nom++)
	{
		fscanf(profile, "%s", mesh_n);

		if (!OpenMesh::IO::read_mesh(noisemeshlist[nom], mesh_n))
		{
			printf("data error");
			return 0;
		}
		halfedgeset[nom].resize(noisemeshlist[nom].n_halfedges());
		for (TriMesh::HalfedgeIter it = noisemeshlist[nom].halfedges_begin(); it != noisemeshlist[nom].halfedges_end(); it++)
		{
			halfedgeset[nom][(*it).idx()].v1 = noisemeshlist[nom].point(noisemeshlist[nom].from_vertex_handle(*it));
			halfedgeset[nom][(*it).idx()].v2 = noisemeshlist[nom].point(noisemeshlist[nom].to_vertex_handle(*it));
		}
		if (noisemeshlist[nom].n_faces() != meshlist[nom].n_faces())
		{
			printf("data error");
			return 0;
		}
	}
	int px[5];  //parameters for gdata 
	// 0,1,2: the index range of output files groups, range(10, 20, 2) = 10, 12, 14, 16, 18
	// 3: the number of LSD in each group
	// 4: the number of LSD in each file
	fscanf(profile, "%s", outputname);  //name and path of training files
	fscanf(profile, "%s", outputflagname); // name and path of ground truth files


	for (int i = 0; i < 5; i++)
		fscanf(profile, "%d", &px[i]);
	if (px[3] % px[4] != 0)
	{

		printf("px3 must be divisible by px4\n");
		return 0;
	}
	int totalcase = 0;
	for (int k0 = px[0]; k0 < px[1]; k0 += px[2])
	{
		totalcase += px[3];
	}
	printf("Output case number: %d\n", totalcase);

	outputcache = new float[px[4] * lsdsize*lsdsize * 3];
	gtcache = new float[px[4] * 3];
	memset(outputcache, 0, px[4] * lsdsize*lsdsize * 3 * sizeof(float));
	memset(gtcache, 0, px[4] * 3*sizeof(float));


	//Make ring, get ground truth normal, noisy normal, face centroid, sigma_s
	for (int nom = 0; nom < numberofmesh; nom++)
	{
		makeRing(meshlist[nom], ringlist[nom], 3);
		getFaceNormal(meshlist[nom], filtered_normals[nom]);
		getFaceNormal(noisemeshlist[nom], noisy_normals[nom]);
		getFaceCentroid(noisemeshlist[nom], face_centroid[nom]);
		sigma_s[nom] = getSigmaS(2, face_centroid[nom], noisemeshlist[nom]) / 8;
		//obtain d_a/p_s
	}


	//flagz=0: the face is far from the boundary
	//flagz=-1: the face is close to the boundary

	for (int nom = 0; nom < numberofmesh; nom++)
	{

		flagz[nom].resize(meshlist[nom].n_faces());
		for (TriMesh::FaceIter v_it = meshlist[nom].faces_begin(); v_it != meshlist[nom].faces_end(); v_it++)
		{
			int index = v_it->idx();
			flagz[nom][index] = 0;
		}
		for (TriMesh::FaceIter v_it = meshlist[nom].faces_begin(); v_it != meshlist[nom].faces_end(); v_it++)
		{
			int index = v_it->idx();
			int count = 0;
			for (TriMesh::FaceFaceIter ff_it = meshlist[nom].ff_begin(TriMesh::FaceHandle(*v_it)); ff_it.is_valid(); ff_it++)
				count++;
			if (count <= 2)
			{
				flagz[nom][index] = -1;
			}
		} //find the faces on the boundary and mark them with -1


		for (TriMesh::FaceIter v_it = meshlist[nom].faces_begin(); v_it != meshlist[nom].faces_end(); v_it++)
		{

			int index = v_it->idx();
			if (flagz[nom][index] <= -1)
				continue;
			for (TriMesh::FaceVertexIter fv_it = meshlist[nom].fv_begin(TriMesh::FaceHandle(index)); fv_it.is_valid(); fv_it++)
			{

				for (TriMesh::VertexFaceIter vf_it = meshlist[nom].vf_begin(*fv_it); vf_it.is_valid(); vf_it++)
				{
					if (flagz[nom][vf_it->idx()] == -1)
					{
						flagz[nom][index] = -2;
						break;
					}
				}
				if (flagz[nom][index] == -2)
					break;
			}

		} //find the faces that their neighbours are on the boundary, and mark them with -2


		for (TriMesh::FaceIter v_it = meshlist[nom].faces_begin(); v_it != meshlist[nom].faces_end(); v_it++)
		{
			int index = v_it->idx();
			if (flagz[nom][index] == -2)
			{
				flagz[nom][index] = -1;
			}
			traindata.push_back(std::pair<int, int>(nom, index));
		}
		// treat -2 and -1 as the same class, collect the mesh index and face index as traindata
	}
	std::random_shuffle(traindata.begin(), traindata.end());
	printf("Total face number: %d\n", traindata.size());

	gsupmet();


	int ttx = -1;

	printf("Generate LSD\n");
	for (int k0 = px[0]; k0 < px[1]; k0 += px[2])
	{

		int count = 0;
		int fcount = 0;
		outputname[strlen(outputname) - 9] = k0 / 10 + '0';
		outputname[strlen(outputname) - 8] = k0 % 10 + '0';
		outputflagname[strlen(outputflagname) - 9] = k0 / 10 + '0';
		outputflagname[strlen(outputflagname) - 8] = k0 % 10 + '0';
		outputname[strlen(outputname) - 6] = '0';
		outputname[strlen(outputname) - 5] = '0';
		outputflagname[strlen(outputflagname) - 6] =  '0';
		outputflagname[strlen(outputflagname) - 5] =  '0';
		if (mt_flag==0)
		{
			for (int k1 = 0; k1 < px[3];)
			{
				ttx++;
				if (ttx == traindata.size())
					ttx = 0;
				int index = traindata[ttx].second;
				int meshidx = traindata[ttx].first;

				if (gLSD(index, noisemeshlist[meshidx], outputcache + count*lsdsize*lsdsize * 3, gtcache + count * 3, sigma_s[meshidx], ringlist[meshidx], filtered_normals[meshidx], halfedgeset[meshidx], noisy_normals[meshidx], face_centroid[meshidx], flagz[meshidx]) == -4)
					continue;
				else
					k1++;


				count++;
				if (count == px[4])
				{
					outputname[strlen(outputname) - 6] = fcount / 10 + '0';
					outputname[strlen(outputname) - 5] = fcount % 10 + '0';
					outputflagname[strlen(outputflagname) - 6] = fcount / 10 + '0';
					outputflagname[strlen(outputflagname) - 5] = fcount % 10 + '0';
					FILE* outfile1 = fopen(outputname, "wb");
					FILE* outfile2 = fopen(outputflagname, "wb");

					fwrite(outputcache, sizeof(float), px[4] * lsdsize*lsdsize * 3, outfile1);
					fwrite(gtcache, sizeof(float), px[4] * 3, outfile2);
					fclose(outfile1);
					fclose(outfile2);
					count = 0;
					fcount++;

					memset(outputcache, 0, px[4] * lsdsize*lsdsize * 3 * sizeof(float));
					memset(gtcache, 0, px[4] * 3 * sizeof(float));

				}
			}
		}
		else
		{
			for (int k1 = 0; k1 < thread_number; k1++)
				thread_p[k1].clear();
			for (int k1 = 0; k1 < px[3];)
			{
				ttx++;
				if (ttx == traindata.size())
					ttx = 0;
				int index = traindata[ttx].second;
				int meshidx = traindata[ttx].first;
				thread_p[count%thread_number].push_back(pid(index, meshidx, count));
				count++;
				k1++;
				if (count == px[4])
				{
					for (int k2 = 0; k2 < thread_number; k2++)
					{
						td[k2] = std::thread(threadprocess, k2);
					}
					for (int k2 = 0; k2 < thread_number; k2++)
					{
						td[k2].join();
					}
					outputname[strlen(outputname) - 6] = fcount / 10 + '0';
					outputname[strlen(outputname) - 5] = fcount % 10 + '0';
					outputflagname[strlen(outputflagname) - 6] = fcount / 10 + '0';
					outputflagname[strlen(outputflagname) - 5] = fcount % 10 + '0';
					FILE* outfile1 = fopen(outputname, "wb");
					FILE* outfile2 = fopen(outputflagname, "wb");

					fwrite(outputcache, sizeof(float), px[4] * lsdsize*lsdsize * 3, outfile1);
					fwrite(gtcache, sizeof(float), px[4] * 3, outfile2);
					fclose(outfile1);
					fclose(outfile2);
					count = 0;
					fcount++;

					memset(outputcache, 0, px[4] * lsdsize*lsdsize * 3 * sizeof(float));
					memset(gtcache, 0, px[4] * 3 * sizeof(float));
					for (int k1 = 0; k1 < thread_number; k1++)
						thread_p[k1].clear();
				}
			}
		}
		if (count > 0)
		{
			printf("error\n");
			return 0;
		}

	}
	delete outputcache;
	delete gtcache;
	return 0;
}

