#include <cmath>
#include "scene.hpp"
#include "vec3.hpp"
#include "ray.hpp"

#define PI 3.141592654f

double uniform_random_number(void){
    return rand()/double(RAND_MAX);
}

scene::scene(int xres, int yres, float fieldOfView, float focalLength, vector3 origin, vector3 lookat, vector3 lookup){
    x_res = xres; //number of pixels in x direction
    y_res= yres; 
    height = 2*(float)focalLength*tan((float)fieldOfView/360.0f *PI/2.0f ); //height in world coords
    width = ((float)xres/(float)yres)*height;
    focal_length = focalLength;
    Camera_origin.setValue(origin.x(), origin.y(), origin.z());
    Camera_lookat.setValue(lookat.x(), lookat.y(), lookat.z());
    Camera_lookup.setValue(lookup.x(), lookup.y(), lookup.z());

    eye_n = vector3::add(Camera_origin, vector3::ScalarMultiply(-1,  Camera_lookat)); //camera coord system
    eye_n.normalize();
    eye_u = vector3::crossproduct(Camera_lookup, eye_n);
    eye_u.normalize();
    eye_v = vector3::crossproduct(eye_n,eye_u); 
    Centre_of_image = vector3::add(Camera_origin,vector3::ScalarMultiply(-focal_length,eye_n));     
    top_left = vector3::add3(Centre_of_image, vector3::ScalarMultiply(-width/2.0f,eye_u), vector3::ScalarMultiply(height/2.0f, eye_v));
    ratio = width/(float)x_res;
}

light::light(float light_length, float light_illumination, vector3 C){ //sets up position of camera
    x_min = C.x()-light_length;
    x_max = C.x()+light_length;
    y_min = C.y()-light_length;
    y_max = C.y()+light_length;
    z_coord = C.z();
    illumination = light_illumination;
    centre.setValue(C.x(), C.y(), C.z()); //set at centre of screen.
    direction.setValue(0.0f,0.0f,-1.0f);
}
vector3 light::PointOnSource(void){ //calculates random point on light source
    float a = uniform_random_number();
    float b = uniform_random_number();
    vector3 tangent_v(0,1,0);
    vector3 tangent_u(1,0,0);
    vector3 Si = vector3::add3(centre, vector3::ScalarMultiply((0.5 - a)*2*x_min,tangent_u), vector3::ScalarMultiply((0.5 -b)*2*x_max,tangent_v));
    return Si;
}
   
sphere_light::sphere_light(vector3 sphere_centre, float sphere_radius){ //spherical light source
    centre.setValue(sphere_centre.x(), sphere_centre.y(), sphere_centre.z());
    radius = sphere_radius;
} 
vector3 sphere_light::PointOnSource(void){
    float a = uniform_random_number()*270.0f; 
    float b = uniform_random_number()*90.0f+90.0f; 
    vector3 Si((float)radius*sin(b/180.0f*PI)*sin(a/180.0f*PI)+centre.x(), (float)radius*sin(b/180.0f*PI)*cos(a/180.0f*PI)+centre.y(),(float)radius*cos(b/180.0f*PI)+centre.z());
    return Si;
}

float sphere_light::intensity(vector3 point){ //inverse square law for light reduction
    vector3 diff = vector3::add(centre, vector3::ScalarMultiply(-1, point));
    float min_diff = (centre.z()-point.z())*(centre.z()-point.z());
    float d = abs(vector3::dotproduct(diff, diff))/min_diff;
    return 10.0f/(d*4.0f*PI);
}