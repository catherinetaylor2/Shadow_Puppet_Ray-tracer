#include <cmath>
#include "scene.hpp"
#include "vec3.hpp"
#include "ray.hpp"

#define PI 3.141592654f

double uniformRandomNumber(void){return rand()/double(RAND_MAX);}

scene::scene(int xres, int yres, float fieldOfView, float focalLength, vector3 origin, vector3 lookat, vector3 lookup){
    xRes = xres; //number of pixels in x direction
    yRes= yres; 
    height = 2*(float)focalLength*tan((float)fieldOfView/360.0f *PI/2.0f ); //height in world coords
    width = ((float)xres/(float)yres)*height;
    focalLength = focalLength;
    cameraOrigin.setValue(origin.x(), origin.y(), origin.z());
    cameraLookat.setValue(lookat.x(), lookat.y(), lookat.z());
    cameraLookup.setValue(lookup.x(), lookup.y(), lookup.z());

    eyeN = vector3::add(cameraOrigin, vector3::ScalarMultiply(-1,  cameraLookat)); //camera coord system
    eyeN.normalize();
    eyeU = vector3::crossproduct(cameraLookup, eyeN);
    eyeU.normalize();
    eyeV = vector3::crossproduct(eyeN,eyeU); 
    centreOfImage = vector3::add(cameraOrigin,vector3::ScalarMultiply(-focalLength,eyeN));     
    topLeft = vector3::add3(centreOfImage, vector3::ScalarMultiply(-width/2.0f,eyeU), vector3::ScalarMultiply(height/2.0f, eyeV));
    ratio = width/(float)xRes;
}

light::light(float lightLength, float lightIllumination, vector3 _centre){ //sets up position of camera
    xMin = _centre.x()-lightLength;
    xMax = _centre.x()+lightLength;
    yMin = _centre.y()-lightLength;
    yMax = _centre.y()+lightLength;
    zCoord = _centre.z();
    illumination = lightIllumination;
    centre.setValue(_centre.x(), _centre.y(), _centre.z()); //set at centre of screen.
    direction.setValue(0.0f,0.0f,-1.0f);
}
vector3 light::PointOnSource(void){ //calculates random point on light source
    float a = uniformRandomNumber();
    float b = uniformRandomNumber();
    vector3 tangentV(0,1,0);
    vector3 tangentU(1,0,0);
    vector3 Si = vector3::add3(centre, vector3::ScalarMultiply((0.5 - a)*2*xMin,tangentU), vector3::ScalarMultiply((0.5 -b)*2*xMax,tangentV));
    return Si;
}