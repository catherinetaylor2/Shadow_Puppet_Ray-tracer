#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"

vector3::vector3(float x, float y, float z){
    x_val=x;
    y_val=y;
    z_val=z;
}
vector3::vector3(){
    x_val=0.0f;
    y_val=0.0f;
    z_val=0.0f;
}
float vector3::get_x(void){
    return x_val;
}
float vector3::get_y(void){
    return y_val;
}
float vector3::get_z(void){
    return z_val;
}
void vector3::setValue(float x, float y, float z){
    x_val=x;
    y_val=y;
    z_val=z;
}
void vector3::normalize(void){
    float sum = sqrt(x_val*x_val + y_val*y_val + z_val*z_val);
    x_val/= sum;
    y_val/= sum;
    z_val/= sum;
}

