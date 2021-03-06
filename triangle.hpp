#ifndef triangle_hpp
#define triangle_hpp

#include "ray.hpp"
#include "vec3.hpp"
#include "binarySearchTree.hpp"
#include "scene.hpp"

#define infinity FLT_MAX;

class triangle{
    public:
        triangle(vector3 V1, vector3 V2, vector3 V3);
        float RayTriangleIntersection(Ray ray);
        static float getPOI(binarySearchTree* root, 
                            float*vertices, 
                            Ray ray,
                            int* faces,
                            int* minValue, //min dist
                            int**k //stores triangle indicies
                            );
        static float getColour( Ray ray, 
                                binarySearchTree* root, 
                                float*vertices, int* faceVertices, int*faceTetxures, float*tetxures, //obj file data
                                unsigned char* puppetTexture, int puppetWidth, int puppetHeight, //texture data
                                vector3 planeNormal, 
                                vector3 LightDir, 
                                float**colours, //stores colours
                                int index //current index
                                );
        static float getColourSoft( Ray rayOuter, Ray rayInner, 
                                    binarySearchTree* root, float*vertices, int* faceVertices, int*faceTetxures, float*tetxures,
                                    unsigned char* puppetTexture, int puppetWidth, int puppetHeight,
                                    vector3 planeNormal,
                                    vector3 LightDir, 
                                    float**colours,
                                    int index,
                                    light innerLight, light outerLight
                                );
        static void getTextureValue(int triangle, //intersection triangle
                                    int* faceVertices, float *vertices, 
                                    Ray ray, 
                                    unsigned char* puppetTexture,
                                    int* faceTetxures, float* tetxures, 
                                    int puppetWidth, int puppetHeight,
                                    float **colour);
    private:
        float planeConstant;
        vector3 vertex1;
        vector3 vertex2;
        vector3 vertex3;
        vector3 normal;
};

#endif