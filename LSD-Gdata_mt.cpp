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
std::vector<double> sigma_s_list;
std::vector<std::vector<ring>> ringlist_list;
std::vector<std::vector<TriMesh::Normal>> filtered_normals_list;
std::vector<std::vector<line>> halfedgeset_list;
std::vector<std::vector<TriMesh::Normal>> noisy_normals_list;
std::vector<std::vector<TriMesh::Point>> face_centroid_list;
std::vector<std::vector<int>> flagz_list;

int gLSD(int index, TriMesh &mesh2, float outputmat[lsdsize*lsdsize*3], float groundtruth[3],
	double sigma_s,
	std::vector<ring> &ringlist,
	std::vector<TriMesh::Normal> &filtered_normals,
	std::vector<line> &halfedgeset,
	std::vector<TriMesh::Normal> &noisy_normals,
	std::vector<TriMesh::Point> &face_centroid,
	std::vector<int> &flagz)
{
	// //obtain n*
	TriMesh::Normal a1 = getAveNormal(ringlist[index], noisy_normals, flagz[index], flagz);

	//obtain polar axis
	TriMesh::Normal startnormal = getPolarAxis(mesh2, index, face_centroid);

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
	int err = samplingNormal(mesh2, index, d2, startnormal, face_centroid, noisy_normals, halfedgeset, sigma_s, lsdsize, outputmat);
	return err;
}

void threadprocess(int p)
{
	for (int i = 0; i < thread_p[p].size(); i++)
	{
		int index = thread_p[p][i].index;
		int meshidx = thread_p[p][i].meshindex;
		int count = thread_p[p][i].count;

		gLSD(index, noisemeshlist[meshidx], outputcache + count*lsdsize*lsdsize * 3, gtcache + count * 3, sigma_s_list[meshidx], ringlist_list[meshidx], filtered_normals_list[meshidx], halfedgeset_list[meshidx], noisy_normals_list[meshidx], face_centroid_list[meshidx], flagz_list[meshidx]);

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
	sigma_s_list.resize(numberofmesh);
	ringlist_list.resize(numberofmesh);
	filtered_normals_list.resize(numberofmesh);
	halfedgeset_list.resize(numberofmesh);
	noisy_normals_list.resize(numberofmesh);
	face_centroid_list.resize(numberofmesh);
	flagz_list.resize(numberofmesh);

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

		ringlist_list[nom].resize(meshlist[nom].n_faces());
		noisy_normals_list[nom].resize(meshlist[nom].n_faces());
		face_centroid_list[nom].resize(meshlist[nom].n_faces());
		filtered_normals_list[nom].resize(meshlist[nom].n_faces());
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
		halfedgeset_list[nom].resize(noisemeshlist[nom].n_halfedges());
		for (TriMesh::HalfedgeIter it = noisemeshlist[nom].halfedges_begin(); it != noisemeshlist[nom].halfedges_end(); it++)
		{
			halfedgeset_list[nom][(*it).idx()].v1 = noisemeshlist[nom].point(noisemeshlist[nom].from_vertex_handle(*it));
			halfedgeset_list[nom][(*it).idx()].v2 = noisemeshlist[nom].point(noisemeshlist[nom].to_vertex_handle(*it));
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
		makeRing(meshlist[nom], ringlist_list[nom], 3);
		getFaceNormal(meshlist[nom], filtered_normals_list[nom]);
		getFaceNormal(noisemeshlist[nom], noisy_normals_list[nom]);
		getFaceCentroid(noisemeshlist[nom], face_centroid_list[nom]);
		sigma_s_list[nom] = getSigmaS(2, face_centroid_list[nom], noisemeshlist[nom]) / 8;
		markBoundaryFaces(meshlist[nom], flagz_list[nom]);
		//obtain d_a/p_s
		for (TriMesh::FaceIter v_it = meshlist[nom].faces_begin(); v_it != meshlist[nom].faces_end(); v_it++){
			int index = v_it -> idx();
			traindata.push_back(std::pair<int, int>(nom, index));
		}
	}
	
	std::random_shuffle(traindata.begin(), traindata.end());
	printf("Total face number: %d\n", traindata.size());

	gsupmat();


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

				if (gLSD(index, noisemeshlist[meshidx], outputcache + count*lsdsize*lsdsize * 3, gtcache + count * 3, sigma_s_list[meshidx], ringlist_list[meshidx], filtered_normals_list[meshidx], halfedgeset_list[meshidx], noisy_normals_list[meshidx], face_centroid_list[meshidx], flagz_list[meshidx]) == -4)
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

