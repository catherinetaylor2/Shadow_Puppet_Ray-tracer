#ifndef ray_hpp
#define ray_hpp

#include <iostream>
#include <vector>
#include "vec3.hpp"

class Ray{
    public:
        Ray(vector3 origin, vector3 direction);
        vector3 getOrigin(void){
            return origin;
        }
        vector3 getDirection(void){
            return direction;
        }
    private:
        vector3 origin;
        vector3 direction;
};
#endif