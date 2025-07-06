#include "LSD.h"
#include <fstream>

std::vector<SampleDirection> local_sample;
std::thread td[thread_number];
float *outputcache;
const int filesize = 10000;
std::vector<line> halfedgeset;
std::vector<TriMesh::Normal> noisy_normals;
std::vector<TriMesh::Point> face_centroid;
std::vector<TriMesh::Normal> filtered_normals;
std::vector<int> flagz;
std::vector<ring> ringlist;
double sigma_s = 0;
const std::string folder = "./virtualData/";

int preprocessing(TriMesh &noisemesh)
{
	ringlist.resize(noisemesh.n_faces());
	noisy_normals.resize(noisemesh.n_faces());
	face_centroid.resize(noisemesh.n_faces());
	filtered_normals.resize(noisemesh.n_faces());
	halfedgeset.resize(noisemesh.n_halfedges());

	// errorflag.resize(noisemesh.n_faces());
	// msave.resize(noisemesh.n_faces());

	for (TriMesh::HalfedgeIter it = noisemesh.halfedges_begin(); it != noisemesh.halfedges_end(); it++)
	{
		halfedgeset[(*it).idx()].v1 = noisemesh.point(noisemesh.from_vertex_handle(*it));
		halfedgeset[(*it).idx()].v2 = noisemesh.point(noisemesh.to_vertex_handle(*it));
	}
	makeRing(noisemesh, ringlist, 3);
	getFaceNormal(noisemesh, noisy_normals);
	getFaceCentroid(noisemesh, face_centroid);
	sigma_s = getSigmaS(2, face_centroid, noisemesh) / 16;

	markBoundaryFaces(noisemesh, flagz);

	return 0;
}

TriMesh noisymesh;
// std::string noisymesh_path = "./strain/bumpy_sphere_n1.obj";
std::string noisymesh_path = "./block.obj";

TriMesh mesh;
std::string mesh_path = "./strain/bumpy_sphere.obj";

std::string num[10] = {"10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};

void savePointsToPly(const std::vector<TriMesh::Point> &points, const std::string &filename)
{
	// Open file for writing
	std::ofstream plyFile(filename);

	if (!plyFile.is_open())
	{
		std::cerr << "Error: Could not open file for writing!" << std::endl;
		return;
	}

	// Write PLY header
	plyFile << "ply\n";
	plyFile << "format ascii 1.0\n";
	plyFile << "element vertex " << points.size() << "\n";
	plyFile << "property float x\n";
	plyFile << "property float y\n";
	plyFile << "property float z\n";
	plyFile << "property uchar red\n";
	plyFile << "property uchar green\n";
	plyFile << "property uchar blue\n";
	plyFile << "element face 0\n"; // We don't have faces in this example, so set to 0
	plyFile << "property list uchar int vertex_indices\n";
	plyFile << "end_header\n";

	// Write vertex data (coordinates and color)
	for (const auto &point : points)
	{
		plyFile << point[0] << " " << point[1] << " " << point[2] << " ";
		plyFile << 255 << " " << 0 << " " << 0 << "\n"; // Set color to red (255, 0, 0)
	}

	plyFile.close();
	std::cout << "PLY file saved to " << filename << std::endl;
}

void saveSelectedFacesToPLY(TriMesh &mesh, const std::vector<int> &traindata, const std::string &filename, int dv)
{
	TriMesh selectedMesh;
	std::string savefile = filename + num[dv / 10 - 1] + ".ply";
	// Add the vertices from the original mesh
	for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		selectedMesh.add_vertex(mesh.point(*v_it));
	}

	// Add the selected faces
	for (int i = 0; i < traindata.size() * dv / 100; ++i)
	{
		int faceIdx = traindata[i];
		TriMesh::FaceHandle fh = mesh.face_handle(faceIdx);

		// Create a new face and add it to the selected mesh
		std::vector<TriMesh::VertexHandle> vertices;
		for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); ++fv_it)
		{
			vertices.push_back(*fv_it);
		}
		selectedMesh.add_face(vertices);
	}

	// Write the selected faces to a PLY file
	if (OpenMesh::IO::write_mesh(selectedMesh, savefile))
	{
		std::cout << "Selected faces saved to " << savefile << std::endl;
	}
	else
	{
		std::cerr << "Error: Could not write mesh to PLY file " << filename << std::endl;
	}
}

int readmesh(TriMesh &mesh, std::string mesh_path)
{

	if (!OpenMesh::IO::read_mesh(mesh, mesh_path))
	{
		std::cerr << "Error: Cannot read mesh from file " << mesh_path << std::endl;
		return 1;
	}
	std::cout << "Mesh Loaded Successfully!" << std::endl;
	return 0;
}

int samplingNormal_test(
	TriMesh &mesh,
	int index,
	const Eigen::Matrix3d &d2,
	const TriMesh::Normal &startnormal,
	const std::vector<TriMesh::Point> &face_centroid,
	const std::vector<TriMesh::Normal> &noisy_normals,
	std::vector<line> &halfedgeset,
	double sigma_s,
	float *outputmat,
	std::vector<int> &sample_face,
	std::vector<MyTraits::Point> &sample_point)
{

	for (int i = 0; i < local_sample.size(); i++)
	{
		// init
		auto &s = local_sample[i];
		double glength = sigma_s * 1.0 * s.radius;
		double clength = 0;
		int centreface = index;

		TriMesh::Point nowpoint = face_centroid[index];
		TriMesh::Normal nownormal = startnormal; // direction of geodesics
		// compute the direction of geodesics
		Eigen::AngleAxisd rotation_vector(
			s.theta, Eigen::Vector3d(
						 noisy_normals[index][0],
						 noisy_normals[index][1],
						 noisy_normals[index][2]));

		Eigen::Vector3d temp3(nownormal[0], nownormal[1], nownormal[2]);
		temp3 = rotation_vector * temp3;
		nownormal = TriMesh::Normal(temp3[0], temp3[1], temp3[2]);
		nownormal.normalize();

		if (i == 0)
		{
			Eigen::Vector3d temp5(
				noisy_normals[centreface][0],
				noisy_normals[centreface][1],
				noisy_normals[centreface][2]);
			temp5 = d2 * temp5;
			outputmat[0] = (float)temp5[0];
			outputmat[1] = (float)temp5[1];
			outputmat[2] = (float)temp5[2];
			sample_face.push_back(centreface);
			sample_point.push_back(nowpoint);
			continue;
		}
		std::unordered_set<int> visitedface;
		TriMesh::FaceHandle nowface(index);
		int endflag = 0;
		int halfedgenum = -1;
		while (endflag == 0)
		{
			visitedface.insert(nowface.idx());
			int edgecount = 0;
			int goflag = 0;
			for (auto it = mesh.fh_begin(nowface); it != mesh.fh_end(nowface); it++)
			{
				int nowhalfedge = it->idx();
				if (nowhalfedge == halfedgenum)
					continue; // 防止跳回来

				edgecount++;
				auto temppoint = nowpoint + nownormal * sigma_s * 100;
				auto tempnormal = nownormal;
				tempnormal.normalize();
				TriMesh::Point nextpoint;
				if (!CalculateLineLineIntersection(halfedgeset[nowhalfedge].v1, halfedgeset[nowhalfedge].v2, nowpoint, temppoint, nextpoint, nownormal))
				{
					continue;
				}
				goflag = 1;
				clength += (nextpoint - nowpoint).length();

				if (clength >= glength)
				{
					Eigen::Vector3d temp5(
						noisy_normals[nowface.idx()][0],
						noisy_normals[nowface.idx()][1],
						noisy_normals[nowface.idx()][2]);
					temp5 = d2 * temp5;
					temp5.normalize();

					outputmat[i * 3 + 0] = (float)temp5[0];
					outputmat[i * 3 + 1] = (float)temp5[1];
					outputmat[i * 3 + 2] = (float)temp5[2];
					sample_face.push_back(nowface.idx());
					sample_point.push_back(nextpoint - tempnormal * (clength - glength));
					endflag = -1;
					break;
				}
				// go to next face
				nowpoint = nextpoint;
				halfedgenum = mesh.opposite_halfedge_handle(*it).idx();
				OpenMesh::FaceHandle nextface = mesh.face_handle(mesh.opposite_halfedge_handle(*it));

				if (nextface.idx() == -1 || visitedface.count(nextface.idx()))
				{
					endflag = -2;
					break;
				}

				Eigen::Matrix3d d = Eigen::Quaterniond::FromTwoVectors(
										Eigen::Vector3d(noisy_normals[nowface.idx()][0],
														noisy_normals[nowface.idx()][1],
														noisy_normals[nowface.idx()][2]),
										Eigen::Vector3d(noisy_normals[nextface.idx()][0],
														noisy_normals[nextface.idx()][1],
														noisy_normals[nextface.idx()][2]))
										.toRotationMatrix();
				Eigen::Vector3d temp2(nownormal[0], nownormal[1], nownormal[2]);
				temp2 = d * temp2;
				nownormal = TriMesh::Normal(temp2[0], temp2[1], temp2[2]);
				nownormal.normalize();

				nowface = nextface;
				break;
			}
			if ((edgecount == 2 && !goflag && halfedgenum != -1) || (edgecount == 3 && !goflag && halfedgenum == -1))
			{
				endflag = -4;
				printf("error: %d %d\n", index, i);
				return endflag;
			}
		}
	}
	return 0;
}

int gLSD_test(int index, float outputmat[(lsd_r_size * lsd_t_size + 1) * 3], std::vector<int> &sample_face, std::vector<MyTraits::Point> &sample_point)
{

	// //obtain n*
	TriMesh::Normal a1 = getAveNormal(ringlist[index], noisy_normals, flagz[index], flagz);

	// obtain polar axis
	TriMesh::Normal startnormal = getPolarAxis(noisymesh, index, face_centroid);

	// //obtain rotation matrix and inverse rotation matrix
	Eigen::Matrix3d d2(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(a1.data()[0],
																		  a1.data()[1],
																		  a1.data()[2]),
														  Eigen::Vector3d(1, 0, 0)));

	Eigen::Matrix3d d2r(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(a1.data()[0],
																									 a1.data()[1],
																									 a1.data()[2])));

	// generate LSD

	int err = samplingNormal_test(noisymesh, index, d2, startnormal, face_centroid, noisy_normals, halfedgeset, sigma_s, outputmat, sample_face, sample_point);
	// return err;
	return err;
}

void saveSampleToPLY(const std::vector<SampleDirection> &local_sample, const std::string &filename)
{
	// 打开PLY文件
	std::ofstream outFile(filename);
	if (!outFile)
	{
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return;
	}

	// 写入PLY头部信息
	outFile << "ply\n";
	outFile << "format ascii 1.0\n";
	outFile << "element vertex " << local_sample.size() << "\n";
	outFile << "property float x\n";
	outFile << "property float y\n";
	outFile << "property float z\n";
	outFile << "property uchar red\n";
	outFile << "property uchar green\n";
	outFile << "property uchar blue\n";
	outFile << "end_header\n";

	// 将local_sample中的采样点写入文件
	for (const auto &sample : local_sample)
	{
		double r = sample.radius;	 // 半径
		double theta = sample.theta; // 角度（以弧度为单位）

		// 计算x和y坐标基于极坐标转换公式
		double x = r * cos(theta); // x = r * cos(theta)
		double y = r * sin(theta); // y = r * sin(theta)
		double z = 0.0f;		   // 默认z坐标为0

		// 设置颜色（例如，显眼的红色）
		unsigned char red = 255;
		unsigned char green = 0;
		unsigned char blue = 0;

		outFile << x << " " << y << " " << z << " " << (int)red << " " << (int)green << " " << (int)blue << "\n";
	}

	outFile.close();
	std::cout << "PLY file " << filename << " has been saved." << std::endl;
}

int main()
{
	outputcache = new float[sampling_size * 3];
	readmesh(noisymesh, noisymesh_path);
	// readmesh(mesh, mesh_path);

	preprocessing(noisymesh);

	int n_faces = noisymesh.n_faces();
	std::vector<int> traindata = globalSampling(noisymesh, flagz, n_faces);

	for (int i = 10; i <= 100; i += 10)
	{
		saveSelectedFacesToPLY(noisymesh, traindata, folder + "GlobalSampling", i);
	}

	std::cout << traindata.size() << std::endl;
	std::vector<int> localSampleResult;
	std::vector<MyTraits::Point> localSamplePoint;
	generateLocalSamplingOrder(local_sample);
	gLSD_test(5178, outputcache, localSampleResult, localSamplePoint);
	for (int i = 10; i <= 100; i += 10)
	{
		saveSelectedFacesToPLY(noisymesh, localSampleResult, folder + "LocalSampling", i);
	}
	savePointsToPly(localSamplePoint, folder + "Point" + ".ply");
	saveSampleToPLY(local_sample, folder + "Order.ply");
	return 0;
}