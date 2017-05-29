#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"
#include "ray.hpp"

Ray::Ray(vector3 o, vector3 d){
    origin.setValue(o.get_x(), o.get_y(), o.get_z());
    direction.setValue(d.get_x(), d.get_y(), d.get_z());
}