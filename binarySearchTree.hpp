#ifndef binarySearchTree_hpp
#define binarySearchTree_hpp

#include <iostream>
#include "vec3.hpp"
#include "ray.hpp"
#include <vector>

class binarySearchTree{
    public:
        binarySearchTree *leftNode;
        binarySearchTree *rightNode;
        float parameters [6];
        int*FacesInNode, NumberOfNodeFaces;
        static void buildTree(float* vertices, int* faces, std::vector<binarySearchTree*> *leafNodes, binarySearchTree**root);
        static void traverseTree(binarySearchTree*root, Ray R, std::vector<int> *output);
        static void findLeafNodes(float* vertices, int*faces, int numberOfFaces, std::vector<binarySearchTree*> *leafNodes);
        static void findParameters(int i, float* vertices, int*faces,  std::vector<float> *parameters, std::vector<float> initialParameters);
        static void deleteTree(binarySearchTree* root);
    private:
};
class BoundingBox{
    public:
        BoundingBox(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);  
        bool rayBoxIntersection(Ray R);    
        float getTmin(void){
            return tmin;
        }
        float getTmax(void){
            return tmax;
        }
    private:
        float parameters [6];
        float tmin;
        float tmax;
};

#endif