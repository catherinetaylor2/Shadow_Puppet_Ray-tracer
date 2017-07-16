#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "triangle.hpp"
#include "vec3.hpp"
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

float triangle::ray_triangle_intersection(Ray R){
    float r = vector3::dotproduct(normal, R.get_direction());
    if (fabs(r) < 0.000000001f){
              return 0;
    }
    float t=(plane_constant - vector3::dotproduct(normal,R.get_origin()) )/r;
    vector3 intersection_point = vector3::vec_add(R.get_origin(), vector3::vec_scal_mult(t,  R.get_direction()));

    if(
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex2, vector3::vec_scal_mult(-1, vertex1)), vector3::vec_add(intersection_point, vector3::vec_scal_mult(-1, vertex1))), normal)>=-0.00001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex3, vector3::vec_scal_mult(-1, vertex2)), vector3::vec_add(intersection_point,vector3::vec_scal_mult(-1, vertex2))), normal)>=-0.00001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::vec_add(vertex1, vector3::vec_scal_mult(-1, vertex3)), vector3::vec_add(intersection_point, vector3::vec_scal_mult(-1, vertex3))), normal)>=-0.00001f))
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
    if(B_root.ray_box_intersection(R)==1){
        search_tree::traverse_tree(root, R, &output);
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
            t_values[z-1]=t;
        }
        for (int z=0; z<(*k)[0]; z++){
            if ((t_values[z]>0)){
                vector3 xyz = vector3::vec_add( R.get_origin(), vector3::vec_scal_mult(t_values[z],R.get_direction()));
                if (t_values[z]<t_min){
                    t_min =t_values[z];
                    *min_value = z;
                }
            }
        }
        delete t_values;
    }
    return t_min; //value of t when intersection occurs
 }

void triangle::get_texture_value(int triangle_value, int* FV, float *V, Ray R, unsigned char* dino_tex, int* FT, float* VT, int dino_width, int dino_height, float **colour)
{
    float denominator;

    
    int c_m1, c_m2, c_m3, n1, n2, n3;
    c_m1 = FV[3*triangle_value] -1, c_m2 = FV[3*triangle_value+1]-1, c_m3 = FV[3*triangle_value+2] -1 ;
    vector3 point1(V[3*c_m1], V[3*c_m1+1], V[3*c_m1+2]);
    vector3 point2(V[3*c_m2], V[3*c_m2+1], V[3*c_m2+2]);
    vector3 point3(V[3*c_m3], V[3*c_m3+1], V[3*c_m3+2]);
    vector3 T = vector3::vec_add(R.get_origin(), vector3::vec_scal_mult(-1, point1));
    vector3 E1 = vector3::vec_add(point2, vector3::vec_scal_mult(-1, point1));
    vector3 E2 = vector3::vec_add(point3, vector3::vec_scal_mult(-1, point1));
    denominator = vector3::dotproduct( vector3::crossproduct(R.get_direction(), E2), E1);

    vector3 M (vector3::dotproduct(vector3::crossproduct(T, E1), E2),vector3::dotproduct(vector3::crossproduct(R.get_direction(), E2), T),vector3::dotproduct( vector3::crossproduct(T, E1), R.get_direction()));
    vector3 tuv = vector3::vec_scal_mult(1.0f/(float)denominator, M);
    float *barycentric = new float[3];
    (barycentric)[0] = 1.0f-(tuv.y()+tuv.z());
    (barycentric)[1] = tuv.y();
    (barycentric)[2] = tuv.z();

    int vt1 = FT[3*triangle_value]-1, vt2 = FT[3*triangle_value+1]-1, vt3 = FT[3*triangle_value+2]-1;
    float vt_1x = VT[2*vt1], vt_1y = VT[2*vt1+1], vt_2x = VT[2*vt2], vt_2y = VT[2*vt2+1], vt_3x = VT[2*vt3], vt_3y = VT[2*vt3+1];

     float u_coord, v_coord, alpha, beta, v12r, v12g, v12b, v34r, v34g, v34b;
    int v1x,v1y, v2x, v4y;
    u_coord = (barycentric[0]*vt_1x +barycentric[1]*vt_2x+barycentric[2]*vt_3x)*dino_width;
    v_coord = (barycentric[0]*vt_1y +barycentric[1]*vt_2y+barycentric[2]*vt_3y)*dino_height;
    v1x = (int)floor(u_coord);
    v1y = (int)ceil(v_coord);
    v2x = (int)ceil(u_coord);
    v4y = (int)floor(v_coord);
    if (v1x<0){
        v1x=0;
    }
    if (v2x<0){
        v2x=0;
    }
    if (v1y<0){
        v1y=0;
    }
    if (v4y<0){
        v4y=0;
    }

    alpha = (float)(u_coord - (v2x - v1x)*v1x)/(float) (v2x - v1x);
    beta = (float)(v_coord - (v1y - v4y)*v4y)/(float) (v1y - v4y);

    if (alpha >1){
        alpha=1;
    }
    if(beta>1){
        beta =1;
    }

    v12r = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x] +  alpha*dino_tex[v1y*dino_width*3 + 3*v2x];
    v12g = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x+1] +alpha*dino_tex[v1y*dino_width*3 + 3*v2x+1];
    v12b = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x+2] +alpha*dino_tex[v1y*dino_width*3 + 3*v2x+2];
    v34g =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x+1] +  alpha*dino_tex[v4y*dino_width*3 + 3*v2x+1];
    v34r =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x] +    alpha*dino_tex[v4y*dino_width*3 + 3*v2x];
    v34b =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x+2] +  alpha*dino_tex[v4y*dino_width*3 + 3*v2x+2];

  (*colour)[0] = (1-beta)*v12r + beta*v34r;
  (*colour)[1] = (1-beta)*v12g + beta*v34g;
  (*colour)[2] = (1-beta)*v12b + beta*v34b;
}