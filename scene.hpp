#ifndef scene_hpp
#define scene_hpp

#include "vec3.hpp"
#include "ray.hpp"

class scene{
    public:
        scene(int xres, int yres, float fieldOfView, float focalLength, vector3 origin, vector3 lookat, vector3 lookup);
        int get_x_res(void){
            return x_res;
        }
        int get_y_res(void){
            return y_res;
        }
        float get_width(void){
            return width;
        }
        float get_height(void){
            return height;
        }
        float get_distance_to_image(void){
            return focal_length;
        }
        vector3 get_centre(void){
            return Centre_of_image;
        }
        vector3 get_corner(void){
            return top_left;
        }
        vector3 get_n(void){
            return eye_n;
        }
        vector3 u(void){
            return eye_u;
        }
        vector3 v(void){
            return eye_v;
        }
        float get_ratio(void){
            return ratio;
        }
    private:
        int x_res;
        int y_res;
        float width;
        float height;
        float focal_length;
        float ratio;
        vector3 Camera_origin;
        vector3 Camera_lookat;
        vector3 Camera_lookup;
        vector3 Centre_of_image;
        vector3 top_left;       
        vector3 eye_n;
        vector3 eye_u;
        vector3 eye_v;
};
class light{ //assuming square light
    public:
        light(float light_length,float light_illumination, vector3 C);
        float get_xmin(void){
            return x_min;
        }
        float get_xmax(void){
            return x_max;
        }
        float get_ymin(void){
            return y_min;
        }
        float get_ymax(void){
            return y_max;
        }
        float get_z(void){
            return z_coord;
        }
        float get_illumination(void){
            return illumination;
        }
        vector3 get_centre(void){
            return centre;
        }
        vector3 get_normal(void){
            return direction;
        }
        vector3 PointOnSource(void);
       
    private:
        vector3 centre;
        vector3 direction;
        float x_min;
        float x_max;
        float y_min;
        float y_max;
        float z_coord;
        float illumination;

};
class sphere_light{ //spherical light source
    public:
        sphere_light(vector3 sphere_centre,float sphere_radius);
        vector3 get_centre(void){
            return centre;
        }
        float get_radius(void){
            return radius;
        }
        vector3 PointOnSource(void);
        float intensity(vector3 point);
    private:
        vector3 centre;
        float radius;
        float illumination;

};
#endif