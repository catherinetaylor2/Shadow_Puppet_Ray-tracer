#ifndef ray_hpp
#define ray_hpp

#include <iostream>
#include <vector>
#include "vec3.hpp"
#include "search_tree.hpp"

class Ray{
    public:
        Ray(vector3 origin, vector3 direction);
        vector3 origin;
        vector3 direction;
        vector3 get_origin(void){
            return origin;
        }
        vector3 get_direction(void){
            return direction;
        }
       
};
#endif