#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"
#include "ray.hpp"

Ray::Ray(vector3 o, vector3 d){
    origin.setValue(o.x(), o.y(), o.z());
    direction.setValue(d.x(), d.y(), d.z());
}