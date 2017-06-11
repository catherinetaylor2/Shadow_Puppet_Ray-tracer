#include <cmath>
#include "scene.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include "triangle.hpp"

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

    eye_n = vector3::vec_add(Camera_origin, vector3::vec_scal_mult(-1,  Camera_lookat));
    eye_n.normalize();
    eye_u = vector3::crossproduct(Camera_lookup, eye_n);
    eye_u.normalize();
    eye_v = vector3::crossproduct(eye_n,eye_u); 
    Centre_of_image = vector3::vec_add(Camera_origin,vector3::vec_scal_mult(-focal_length,eye_n));     
    top_left = vector3::vec_add3(Centre_of_image, vector3::vec_scal_mult(-width/2.0f,eye_u), vector3::vec_scal_mult(height/2.0f, eye_v));
    ratio = width/(float)x_res;
}

light::light(float light_length,  float z, float light_illumination){
    x_min = -light_length;
    x_max = light_length;
    y_min = -light_length;
    y_max = light_length;
    z_coord = z;
    illumination = light_illumination;
    centre.setValue(0.0f,0.0f, z_coord);
    direction.setValue(0.0f,0.0f,1.0f);
}
bool light::ray_intersection(Ray R, triangle upper, triangle lower){
    bool t_1 = upper.ray_triangle_intersection(R);
    bool t_2 = lower.ray_triangle_intersection(R);
    if((t_1==1)||(t_2)==1){
        return 1;
    }
    else{
        return 0;
    }
}
float light::DiffuseValue( vector3 normal, vector3 light_direction){
    if (vector3::dotproduct(normal,light_direction)>0){        
        return vector3::dotproduct(normal,light_direction);
    }
    else{
        return 0;
    }
}
   
sphere_light::sphere_light(vector3 sphere_centre, float sphere_radius){
    centre.setValue(sphere_centre.x(), sphere_centre.y(), sphere_centre.z());
    radius = sphere_radius;
} 