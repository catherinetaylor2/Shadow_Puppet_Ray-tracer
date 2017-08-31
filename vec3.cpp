#include <cstdlib>
#include <iostream>
#include <cmath>
#include "vec3.hpp"

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

