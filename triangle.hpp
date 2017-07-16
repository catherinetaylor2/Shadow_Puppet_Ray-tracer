#ifndef triangle_hpp
#define triangle_hpp

#include "ray.hpp"
#include "vec3.hpp"
#include "search_tree.hpp"
#include "scene.hpp"

#define infinity FLT_MAX;

class triangle{
    public:
        triangle(vector3 V1, vector3 V2, vector3 V3);
        vector3 get_triangle_normal(void){
            return normal;
        }
        float ray_triangle_intersection(Ray R);
        static float intersection_point(search_tree* root, float*vertices, Ray R, int* faces, int* min_value, int**k);
        static float intersection_value(Ray R, search_tree* root, float*vertices, int* FV, int*FT, float*VT, unsigned char* dino_tex, int dino_width, int dino_height, vector3 plane_n, vector3 L, float**colours, int index);
        static float intersection_value_s(Ray Rl, Ray R, search_tree* root, float*vertices, int* FV, int*FT, float*VT, unsigned char* dino_tex, int dino_width, int dino_height, vector3 plane_n, vector3 L, float**colours, int index, light inner_light, light outer_light);
        static void get_texture_value(int triangle, int* FV, float *V, Ray R, unsigned char* dino_tex, int* FT, float* VT, int dino_width, int dino_height, float **colour);
    private:
        float plane_constant;
        vector3 vertex1;
        vector3 vertex2;
        vector3 vertex3;
        vector3 normal;
};

#endif