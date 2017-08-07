#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"

vector3::vector3(float x, float y, float z){
    xVal=x;
    yVal=y;
    zVal=z;
}
vector3::vector3(){
    xVal=0.0f;
    yVal=0.0f;
    zVal=0.0f;
}
float vector3::x(void){
    return xVal;
}
float vector3::y(void){
    return yVal;
}
float vector3::z(void){
    return zVal;
}
void vector3::setValue(float x, float y, float z){
    xVal=x;
    yVal=y;
    zVal=z;
}
void vector3::normalize(void){
    float sum = sqrt(xVal*xVal + yVal*yVal + zVal*zVal);
    xVal/= sum;
    yVal/= sum;
    zVal/= sum;
}

