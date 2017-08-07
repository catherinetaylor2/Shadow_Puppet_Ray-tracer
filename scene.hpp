#ifndef scene_hpp
#define scene_hpp

#include "vec3.hpp"
#include "ray.hpp"

class scene{
    public:
        scene(int xres, int yres, float fieldOfView, float focalLength, vector3 origin, vector3 lookat, vector3 lookup);
        int getXRes(void) {return xRes;}
        int getYRes(void) {return yRes;}
        float getWidth(void) {return width;}
        float getHeight(void) {return height;}
        float getDistanceToImage(void) {return focalLength;}
        vector3 getCentre(void) {return centreOfImage;}
        vector3 getCorner(void) {return topLeft;}
        vector3 n(void) {return eyeN;}
        vector3 u(void) {return eyeU;}
        vector3 v(void) {return eyeV;}
        float getRatio(void) {return ratio;}
    private:
        int xRes, yRes;
        float width, height;
        float focalLength;
        float ratio;
        vector3 cameraOrigin, cameraLookat, cameraLookup;
        vector3 centreOfImage;
        vector3 topLeft;       
        vector3 eyeN, eyeU, eyeV;
};
class light{ //assuming square light
    public:
        light(float lightLength,float lightIllumination, vector3 _centre);
        float getXmin(void) {return xMin;}
        float get_xmax(void) {return xMax;}
        float getYmin(void) {return yMin;}
        float getYmax(void) {return yMax;}
        float z(void) {return zCoord;}
        float get_illumination(void) {return illumination;}
        vector3 getCentre(void) {return centre;}
        vector3 get_normal(void) {return direction;}
        vector3 PointOnSource(void);
    private:
        vector3 centre;
        vector3 direction;
        float xMin, xMax;
        float yMin, yMax;
        float zCoord;
        float illumination;
};
#endif