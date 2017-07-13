#ifndef triangle_hpp
#define triangle_hpp

#include "ray.hpp"
#include "vec3.hpp"
#include "search_tree.hpp"

#define infinity FLT_MAX;

class triangle{
    public:
        triangle(vector3 V1, vector3 V2, vector3 V3);
        vector3 get_triangle_normal(void){
            return normal;
        }
        float ray_triangle_intersection(Ray R);
        static float intersection_point(search_tree* root, float*vertices, Ray R, int* faces, int* min_value, int**k);
    private:
        float plane_constant;
        vector3 vertex1;
        vector3 vertex2;
        vector3 vertex3;
        vector3 normal;
};

#endif