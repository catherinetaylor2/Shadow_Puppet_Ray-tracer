#include <cmath>
#include "scene.hpp"
#include "vec3.hpp"
#include "ray.hpp"

#define PI 3.141592654f

scene::scene(int xres, int yres, float fieldOfView, float focalLength, vector3 origin, vector3 lookat, vector3 lookup){
    x_res = xres;
    y_res= yres;
    height = 2*(float)focalLength*tan((float)fieldOfView/360.0f *PI/2.0f );
    width = ((float)xres/(float)yres)*height;
    focal_length = focalLength;
    Camera_origin.setValue(origin.x(), origin.y(), origin.z());
    Camera_lookat.setValue(lookat.x(), lookat.y(), lookat.z());
    Camera_lookup.setValue(lookup.x(), lookup.y(), lookup.z());

    vector3 w = vector3::vec_add(Camera_origin, vector3::vec_scal_mult(-1,  Camera_lookat));
    w.normalize();
    vector3 u = vector3::crossproduct(Camera_lookup, w);
    u.normalize();
    vector3 v = vector3::crossproduct(w,u); 
    Centre_of_image = vector3::vec_add(Camera_origin,vector3::vec_scal_mult(-focal_length,w));     
    top_left = vector3::vec_add3(Centre_of_image, vector3::vec_scal_mult(-width/2.0f,u), vector3::vec_scal_mult(height/2.0f,v));
}

light::light(vector3 light_centre, float light_radius, float light_illumination){
    centre.setValue(light_centre.x(), light_centre.y(), light_centre.z());
    radius = light_radius;
    illumination = light_illumination;
}