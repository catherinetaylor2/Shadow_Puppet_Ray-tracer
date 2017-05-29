#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "triangle.hpp"
#include "vec3.hpp"
#include "scene.hpp"
#include "ray.hpp"
#include "search_tree.hpp"

triangle::triangle(vector3 V1, vector3 V2, vector3 V3){
    vertex1.setValue(V1.x(), V1.y(), V1.z());
    vertex2.setValue(V2.x(), V2.y(), V2.z());
    vertex3.setValue(V3.x(), V3.y(), V3.z());

    vector3 N = vector3::crossproduct(vector3::vec_add(vertex2, vector3::vec_scal_mult(-1, vertex1)), vector3::vec_add(vertex3, vector3::vec_scal_mult(-1, vertex1)));
    N.normalize();
    normal.setValue(N.x(), N.y(), N.z()); 
    plane_constant = vector3::dotproduct(normal, vertex1);
}

bool triangle::ray_triangle_intersection(Ray R){
    float r = vector3::dotproduct(normal, R.get_direction());
    if (fabs(r) < 0.000000001f){
              return 0;
    }
    float t=(plane_constant - vector3::dotproduct(normal,R.get_origin()) )/r;
    vector3 intersection_point = vector3::vec_add(R.get_origin(), vector3::vec_scal_mult(t,  R.get_direction()));

    if(
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex2, vector3::vec_scal_mult(-1, vertex1)), vector3::vec_add(intersection_point, vector3::vec_scal_mult(-1, vertex1))), normal)>=-0.0000001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex3, vector3::vec_scal_mult(-1, vertex2)), vector3::vec_add(intersection_point,vector3::vec_scal_mult(-1, vertex2))), normal)>=-0.0000001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex1, vector3::vec_scal_mult(-1, vertex3)), vector3::vec_add(intersection_point, vector3::vec_scal_mult(-1, vertex3))), normal)>=-0.0000001f))
    {
    return t;
    }
    return 0;
}
 float triangle::intersection_point(search_tree* root, float*vertices, Ray R, int* faces, int* min_value, int**k){
    Bounding_box B_root(root->parameters[0],root->parameters[1], root->parameters[2],root->parameters[3],root->parameters[4],root->parameters[5]);
    std::vector<int> output;
    float t_min = infinity;
	int index1, index2, index3;
    output.clear();
    if(B_root.ray_box_intersection(R.get_origin(), R.get_direction())==1){
        search_tree::traverse_tree(root, R.get_origin(), R.get_direction(), &output);
    }
    *k = new int[output.size()+1];
    (*k)[0] = -1;
    if (output.size()>1){
        (*k)[0]=output.size();
        for(int g=1; g<(*k)[0]+1;g++){
            (*k)[g] = output[g-1];
        }
    }
    if( ((*k)[0]!=-1)&&((*k)[0]>0)){
        float* t_values = new float[(*k)[0]], t;
        int index;
        for (int z=1; z<(*k)[0]+1; z++){
            index = (*k)[z];
            index1 = faces[3*index] -1, index2 = faces[3*index+1]-1, index3 = faces[3*index+2] -1 ;
            vector3 V1(vertices[3*index1], vertices[3*index1+1], vertices[3*index1+2]);
            vector3 V2(vertices[3*index2], vertices[3*index2+1], vertices[3*index2+2]);
            vector3 V3(vertices[3*index3], vertices[3*index3+1], vertices[3*index3+2]);
            triangle tri(V1, V2, V3);
            t = tri.ray_triangle_intersection(R);
            t_values[z-1]=t;
        }
        for (int z=0; z<(*k)[0]; z++){
            if ((t_values[z]>0)){
                vector3 xyz = vector3::vec_add( R.get_origin(), vector3::vec_scal_mult(t_values[z],R.get_direction()));
                if (xyz.z()<t_min){
                    t_min = xyz.z();
                    *min_value = z;
                }
            }
        }
        delete t_values;
    }
    return t_min;
 }