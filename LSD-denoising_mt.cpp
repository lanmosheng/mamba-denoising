#include"LSD.h"

int supmat[lsdsize][lsdsize][3];
std::vector<SampleDirection> gloabl_sample;
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
const int filesize = 10000;

std::vector<pid> thread_p[thread_number];

double sigma_s = 0;
std::vector<ring> ringlist;
TriMesh noisemesh;
std::vector<line> halfedgeset;
std::vector<TriMesh::Normal> noisy_normals;
std::vector<TriMesh::Point> face_centroid;
std::vector<TriMesh::Normal> filtered_normals;
std::vector<int> flagz;

std::vector<Eigen::Matrix3d> msave;
std::vector<int> errorflag;
std::vector<FILE*> filepo;

int gLSD(int index, float outputmat[lsdsize*lsdsize * 3])
{

	// //obtain n*
	TriMesh::Normal a1 = getAveNormal(ringlist[index], noisy_normals, flagz[index], flagz);

	//obtain polar axis
	TriMesh::Normal startnormal = getPolarAxis(noisemesh, index, face_centroid);

	//obtain rotation matrix and inverse rotation matrix             
	Eigen::Matrix3d d2(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2]), Eigen::Vector3d(1, 0, 0)));

	Eigen::Matrix3d d2r(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2])));

	msave[index] = d2r;

	//generate LSD

	int err = samplingNormal(noisemesh, index, d2, startnormal, face_centroid, noisy_normals, halfedgeset, sigma_s, lsdsize, outputmat);
	return err;
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
	
	markBoundaryFaces(noisemesh, flagz);

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

	gsupmat();

	//read noisy meshes 
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
			std::vector<int> sorted_face_order = globalSampling(noisemesh, flagz, n_faces);
			
			for(int i = 0; i < n_faces; i++){
				int face_idx = sorted_face_order[i];
				thread_p[count % thread_number].push_back(pid(face_idx, count));
				count++;

				if(count == filesize || i == n_faces - 1){
					for(int k2 = 0; k2 < thread_number; k2++){
						td[k2] = std::thread(threadprocess, k2);
					}
					for(int k2 = 0; k2 < thread_number; k2++){
						td[k2].join();
					}
					fwrite(outputcache, sizeof(float), count * lsdsize * lsdsize * 3, filepo[fcount]);

					count = 0;
					fcount++;
					memset(outputcache, 0, filesize * lsdsize * lsdsize * 3 * sizeof(float));
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

