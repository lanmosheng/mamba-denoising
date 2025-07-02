#include"LSD.h"

#include <GLFW/glfw3.h>

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

void saveSelectedFacesToPLY(TriMesh& mesh, const std::vector<int>& traindata, std::string& filename, int dv) {
    TriMesh selectedMesh;
    filename = filename + num[dv/10-1] + ".ply";
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
    if (OpenMesh::IO::write_mesh(selectedMesh, filename)) {
        std::cout << "Selected faces saved to " << filename << std::endl;
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
int main(){
    outputcache = new float[filesize * sampling_size * 3];
    readmesh(noisymesh, noisymesh_path);
    // readmesh(mesh, mesh_path);
 
    preprocessing(noisymesh);

    int n_faces = noisymesh.n_faces();
    std::vector<int> traindata = globalSampling(noisymesh, flagz, n_faces);
    
    std::string filen = "sample";
    for(int i=10;i<=100;i+=10){
        saveSelectedFacesToPLY(noisymesh, traindata, filen, i);
        filen = "sample";
    }
    return 0;
}