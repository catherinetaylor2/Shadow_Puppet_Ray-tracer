#include <iostream>
#include "binarySearchTree.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include <vector>
#include <algorithm>

#define infinity FLT_MAX

 BoundingBox::BoundingBox(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){ //define the box with min/max values
     parameters[0] = xmin;
     parameters[1] = xmax;
     parameters[2] = ymin;
     parameters[3] = ymax;
     parameters[4] = zmin;
     parameters[5] = zmax;
 }

bool BoundingBox::rayBoxIntersection(Ray R){ //find boxes which ray intersects with
    float tmin_y, tmax_y, tmin_z, tmax_z;
    vector3 invDirection(1/(R.getDirection().x()), 1/R.getDirection().y(), 1/R.getDirection().z());
    tmin = (parameters[0]- R.getOrigin().x())*invDirection.x();
    tmax = (parameters[1]-R.getOrigin().x())*invDirection.x();
    tmin_y = (parameters[2]- R.getOrigin().y())*invDirection.y();
    tmax_y = (parameters[3]- R.getOrigin().y())*invDirection.y();

     if (tmin > tmax)std::swap(tmin, tmax); 
     if (tmin_y > tmax_y)std::swap(tmin_y, tmax_y);  
     if ((tmin > tmax_y)||(tmin_y>tmax)){
         return 0;
     }
     if (tmin_y > tmin){
        tmin = tmin_y;
     }
     if (tmax_y < tmax){ //only test in x and y as puppet is 2D
         tmax = tmax_y;
     }        
     return 1;     
 }

void binarySearchTree::traverseTree(binarySearchTree*root, Ray R, std::vector<int> *output){ //traverse tree postorder
    if(((root->leftNode!=nullptr))||((root->rightNode!=nullptr))){ //if not null go left, ow go right
        BoundingBox rightBox(root->rightNode->parameters[0],root->rightNode->parameters[1], root->rightNode->parameters[2],root->rightNode->parameters[3],
                             root->rightNode->parameters[4],root->rightNode->parameters[5]);
        BoundingBox leftBox(root->leftNode->parameters[0],root->leftNode->parameters[1], root->leftNode->parameters[2],root->leftNode->parameters[3],
                            root->leftNode->parameters[4],root->leftNode->parameters[5]);
        if ((root->rightNode!=nullptr)&&(rightBox.rayBoxIntersection(R)==1)) {
            traverseTree(root->rightNode,R, output);
        }  
        if ((root->leftNode!=nullptr)&&(leftBox.rayBoxIntersection(R)==1)){
            traverseTree(root->leftNode, R, output);
        }   
    } 
    if(((root->leftNode==nullptr))&&((root->rightNode==nullptr))){ //if null visit root
        BoundingBox B_root(root->parameters[0],root->parameters[1], root->parameters[2],root->parameters[3],root->parameters[4],root->parameters[5]);
        if((B_root.rayBoxIntersection(R)==1)){ //if intersection then save
            for (int i = 0; i<root->NumberOfNodeFaces; i++){
                (*output).push_back( root->FacesInNode[i]);
            }         
        }   
    }
}

void binarySearchTree::findParameters(int i, float* vertices,int*faceVertices, std::vector<float> *parameters, std::vector<float> initialParameters){ //use vertex values to find max and min values
    float xmin = initialParameters[0], ymin = initialParameters[2], zmin = initialParameters[4],
          xmax=initialParameters[1], ymax=initialParameters[3], zmax = initialParameters[5];
    for(int j=0; j<3; j++){
        if (vertices[3*faceVertices[3*i+j]]< xmin){
        xmin = vertices[3*faceVertices[3*i+j]];
        }
        if (vertices[3*faceVertices[3*i+j]+1]< ymin){
        ymin = vertices[3*faceVertices[3*i+j]+1];
        }
        if (vertices[3*faceVertices[3*i+j]+2]< zmin){
        zmin = vertices[3*faceVertices[3*i+j]+2];
        }
        if (vertices[3*faceVertices[3*i+j]]> xmax){
        xmax = vertices[3*faceVertices[3*i+j]];
        }
        if (vertices[3*faceVertices[3*i+j]+1]> ymax){
        ymax =vertices[3*faceVertices[3*i+j]+1];
        }
        if (vertices[3*faceVertices[3*i+j]+2]> zmax){
        zmax = vertices[3*faceVertices[3*i+j]+2];
        }
    }   
    *parameters = {xmin, xmax, ymin, ymax, zmin, zmax};
}

void binarySearchTree::findLeafNodes(float* vertices, int*faceVertices, int number_of_faceVertices, std::vector<binarySearchTree*> *leafNodes){ 
    //create vector full of leaf nodes
    std::vector<float> parameters;
    std::vector<float> initialParameters = {infinity, -1*infinity, infinity, -1*infinity, infinity, -1*infinity};
    for(int i = 0; i<number_of_faceVertices; i++){ 
        binarySearchTree* leaf = new binarySearchTree;
        leaf->NumberOfNodeFaces = 1; 
        leaf->FacesInNode= new int [1];
        leaf->FacesInNode[0]=i;
        parameters.clear();
        binarySearchTree::findParameters(i, vertices, faceVertices, &parameters, initialParameters);
        for(int j =0; j<6;j++){
            leaf->parameters[j] = parameters[j];
        }
        leaf->leftNode=nullptr;
        leaf->rightNode=nullptr;
        (*leafNodes).push_back(leaf);
    }    
}

void binarySearchTree::buildTree(float* vertices, int* faceVertices, std::vector<binarySearchTree*>* leafNodes, binarySearchTree**root){ 
    //use leaf node vector to build bst.
    std::vector<float> parameters;
   
    while((*leafNodes).size()>1){ //while vector has more than one element
        parameters.clear();
        binarySearchTree* temp = new binarySearchTree; //create temp root node
        binarySearchTree* first = (*leafNodes)[0]; //pop first two elements
        binarySearchTree* second = (*leafNodes)[1];
        temp->leftNode =(*leafNodes)[0]; //set these trees as left and right values of temp
        temp->rightNode = (*leafNodes)[1];
        temp->NumberOfNodeFaces = first->NumberOfNodeFaces+second->NumberOfNodeFaces;
        (*leafNodes).erase ((*leafNodes).begin(),(*leafNodes).begin()+2);
        temp->FacesInNode = new int [temp->NumberOfNodeFaces];
        for(int i=0; i<first->NumberOfNodeFaces; i++){
            temp->FacesInNode[i] = first->FacesInNode[i];
        }
        for(int i = first->NumberOfNodeFaces; i<first->NumberOfNodeFaces+second->NumberOfNodeFaces;i++){
            temp->FacesInNode[i] = second->FacesInNode[i-first->NumberOfNodeFaces];
        }
        std::vector<float> initialParameters = {infinity, -1*infinity, infinity, -1*infinity, infinity, -1*infinity};
        for(int i =0; i<temp->NumberOfNodeFaces; i++){
            binarySearchTree::findParameters(temp->FacesInNode[i], vertices, faceVertices, &parameters, initialParameters);
            initialParameters = {parameters[0],parameters[1],parameters[2], parameters[3], parameters[4], parameters[5]};              
        }
        for(int j =0; j<6;j++){
             temp->parameters[j] = parameters[j];
        }
        (*leafNodes).push_back(temp); //add temp value to end of vector
    }
    *root = (*leafNodes)[0];
    return;
}

void binarySearchTree::deleteTree(binarySearchTree* root){ //delete all nodes in binary serach tree;
    if((root->leftNode)!=nullptr){
        deleteTree((root)->leftNode);
    }
     if((root->rightNode)!=nullptr){
        deleteTree((root)->rightNode);
     }
    delete root->FacesInNode;
    delete root;
    root = nullptr;   
}