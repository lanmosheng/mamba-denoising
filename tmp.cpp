#include"LSD.h"
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