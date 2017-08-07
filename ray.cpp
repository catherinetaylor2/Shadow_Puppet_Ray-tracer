#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"
#include "ray.hpp"

#define infinity FLT_MAX;

Ray::Ray(vector3 _origin, vector3 _direction){
    origin.setValue(_origin.x(), _origin.y(), _origin.z());
    direction.setValue(_direction.x(), _direction.y(), _direction.z());
}