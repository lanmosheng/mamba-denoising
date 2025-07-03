#include"LSD.h"


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

int preprocessing(TriMesh& noisemesh)
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
	sigma_s = getSigmaS(2, face_centroid, noisemesh) / 8;
	
	markBoundaryFaces(noisemesh, flagz);

	return 0;
}


TriMesh noisymesh;
std::string noisymesh_path = "./strain/bumpy_sphere_n1.obj";
TriMesh mesh;
std::string mesh_path = "./strain/bumpy_sphere.obj";

std::string num[10]={"10","20","30","40","50","60","70","80","90","100"};

void saveSelectedFacesToPLY(TriMesh& mesh, const std::vector<int>& traindata, const std::string& filename, int dv) {
    TriMesh selectedMesh;
    std::string savefile = filename + num[dv/10-1] + ".ply";
    // Add the vertices from the original mesh
    for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        selectedMesh.add_vertex(mesh.point(*v_it));
    }

    // Add the selected faces
    for (int i = 0; i < traindata.size()*dv/100; ++i) {
        int faceIdx = traindata[i];
        TriMesh::FaceHandle fh = mesh.face_handle(faceIdx);

        // Create a new face and add it to the selected mesh
        std::vector<TriMesh::VertexHandle> vertices;
        for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh); ++fv_it) {
            vertices.push_back(*fv_it);
        }
        selectedMesh.add_face(vertices);
    }

    // Write the selected faces to a PLY file
    if (OpenMesh::IO::write_mesh(selectedMesh, savefile)) {
        std::cout << "Selected faces saved to " << savefile << std::endl;
    } else {
        std::cerr << "Error: Could not write mesh to PLY file " << filename << std::endl;
    }
}

int readmesh(TriMesh& mesh, std::string mesh_path){
    
    if(!OpenMesh::IO::read_mesh(mesh,mesh_path)){
        std::cerr << "Error: Cannot read mesh from file " << mesh_path << std::endl;
        return 1;
    }
    std::cout<<"Mesh Loaded Successfully!"<<std::endl;
    return 0;
}


int samplingNormal_test(
	TriMesh& mesh,
    int index,
	const Eigen::Matrix3d& d2,
    const TriMesh::Normal& startnormal,
    const std::vector<TriMesh::Point>& face_centroid,
    const std::vector<TriMesh::Normal>& noisy_normals,
    std::vector<line>& halfedgeset,
    double sigma_s,
    float* outputmat,
    std::vector<int> &ret)
{
    std::cout<<"Local Sample size:" << local_sample.size()<<std::endl;
	
	for (int i = 0; i < local_sample.size(); i++){
		//init
		auto &s = local_sample[i];
		double glength = sigma_s * 1.0 * s.radius;
		double clength = 0;
		int centreface = index;

		TriMesh::Point nowpoint = face_centroid[index];
		TriMesh::Normal nownormal = startnormal;//direction of geodesics 
		//compute the direction of geodesics
		Eigen::AngleAxisd rotation_vector(
			s.theta, Eigen::Vector3d(
				noisy_normals[index][0],
				noisy_normals[index][1],
				noisy_normals[index][2]
			)
		);

		Eigen::Vector3d temp3(nownormal[0],nownormal[1],nownormal[2]);
		temp3 = rotation_vector * temp3;
		nownormal = TriMesh::Normal(temp3[0],temp3[1],temp3[2]);
		nownormal.normalize();
		
		if(i == 0){
			Eigen::Vector3d temp5(
				noisy_normals[centreface][0],
				noisy_normals[centreface][1],
				noisy_normals[centreface][2]	
			);
			temp5 = d2 * temp5;
			outputmat[0] = (float)temp5[0];
			outputmat[1] = (float)temp5[1];
			outputmat[2] = (float)temp5[2];
            ret.push_back(centreface);
			continue;
		}
		std::unordered_set<int> visitedface;
		TriMesh::FaceHandle nowface(index);
		int endflag = 0;
		int halfedgenum = -1;
		while(endflag == 0){
			visitedface.insert(nowface.idx());
			int edgecount = 0;
			int goflag = 0;
			for(auto it = mesh.fh_begin(nowface); it != mesh.fh_end(nowface); it++){
				int nowhalfedge = it->idx();
				if(nowhalfedge == halfedgenum) continue; //防止跳回来

				edgecount++;
				auto temppoint = nowpoint + nownormal * sigma_s * 100;
				TriMesh::Point nextpoint;
				if(!CalculateLineLineIntersection(halfedgeset[nowhalfedge].v1,halfedgeset[nowhalfedge].v2,nowpoint, temppoint,nextpoint, nownormal)){
					continue;
				}
				goflag = 1;
				clength += (nextpoint - nowpoint).length();

				if(clength >= glength){
					Eigen::Vector3d temp5(
						noisy_normals[nowface.idx()][0],
						noisy_normals[nowface.idx()][1],
						noisy_normals[nowface.idx()][2]
					);
					temp5 = d2 * temp5;
					temp5.normalize();
	
					outputmat[i * 3 + 0] = (float)temp5[0];
					outputmat[i * 3 + 1] = (float)temp5[1];
					outputmat[i * 3 + 2] = (float)temp5[2];
                    ret.push_back(nowface.idx());
					endflag = -1;
					break;
				}
				//go to next face
				nowpoint = nextpoint;
				halfedgenum = mesh.opposite_halfedge_handle(*it).idx();
				OpenMesh::FaceHandle nextface = mesh.face_handle(mesh.opposite_halfedge_handle(*it));

				if(nextface.idx() == -1 || visitedface.count(nextface.idx())){
					endflag = -2;
					break;
				}

				Eigen::Matrix3d d = Eigen::Quaterniond::FromTwoVectors(
					Eigen::Vector3d(noisy_normals[nowface.idx()][0],
									noisy_normals[nowface.idx()][1],
									noisy_normals[nowface.idx()][2]),
					Eigen::Vector3d(noisy_normals[nextface.idx()][0],
									noisy_normals[nextface.idx()][1],
									noisy_normals[nextface.idx()][2])
				).toRotationMatrix();
				Eigen::Vector3d temp2(nownormal[0], nownormal[1], nownormal[2]);
				temp2 = d * temp2;
				nownormal = TriMesh::Normal(temp2[0], temp2[1], temp2[2]);
				nownormal.normalize();

				nowface = nextface;
				break;
			}
			if ((edgecount == 2 && !goflag && halfedgenum != -1) ||(edgecount == 3 && !goflag && halfedgenum == -1)) {
				endflag = -4;
				printf("error: %d %d\n", index, i);
				return endflag;
			}
		}
	}
	
	return 0;
}

int gLSD_test(int index, float outputmat[(lsd_r_size * lsd_t_size + 1) * 3], std::vector<int> &ret)
{

	// //obtain n*
	TriMesh::Normal a1 = getAveNormal(ringlist[index], noisy_normals, flagz[index], flagz);

	//obtain polar axis
	TriMesh::Normal startnormal = getPolarAxis(noisymesh, index, face_centroid);

	// //obtain rotation matrix and inverse rotation matrix             
	Eigen::Matrix3d d2(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2]), Eigen::Vector3d(1, 0, 0)));

	Eigen::Matrix3d d2r(Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(a1.data()[0],
		a1.data()[1],
		a1.data()[2])));

	//generate LSD

	int err = samplingNormal_test(noisymesh, index, d2, startnormal, face_centroid, noisy_normals, halfedgeset, sigma_s, outputmat, ret);
	//return err;
	return err;
}

int main(){
    outputcache = new float[sampling_size * 3];
    readmesh(noisymesh, noisymesh_path);
    // readmesh(mesh, mesh_path);
    
    preprocessing(noisymesh);
    
    int n_faces = noisymesh.n_faces();
    std::vector<int> traindata = globalSampling(noisymesh, flagz, n_faces);
    
    for(int i=10;i<=100;i+=10){
            saveSelectedFacesToPLY(noisymesh, traindata, folder + "GlobalSampling", i);
    }
        
        
    std::cout<<traindata.size()<<std::endl;
    std::vector<int> localSampleResult;
    std::cout<< localSampleResult.size() << std::endl;
    generateLocalSamplingOrder(local_sample);
	gLSD_test(traindata[0], outputcache, localSampleResult);
	for(int i=10;i<=100;i+=10){
		saveSelectedFacesToPLY(noisymesh, localSampleResult, folder + "LocalSampling", i);
	}
    return 0;
}