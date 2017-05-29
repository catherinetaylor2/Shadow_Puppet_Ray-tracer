#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "triangle.hpp"
#include "vec3.hpp"
#include "scene.hpp"
#include "ray.hpp"

triangle::triangle(vector3 V1, vector3 V2, vector3 V3){
    vertex1.setValue(V1.x(), V1.y(), V1.z());
    vertex2.setValue(V2.x(), V2.y(), V2.z());
    vertex3.setValue(V3.x(), V3.y(), V3.z());

    vector3 N = vector3::crossproduct(vector3::vec_add(vertex2, vector3::vec_scal_mult(-1, vertex1)), vector3::vec_add(vertex3, vector3::vec_scal_mult(-1, vertex1)));
    N.normalize();
    normal.setValue(N.x(), N.y(), N.z()); 
    plane_constant = vector3::dotproduct(normal, vertex1);
}