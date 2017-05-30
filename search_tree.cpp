#include <iostream>
#include "search_tree.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include <vector>
#include <algorithm>

#define infinity FLT_MAX

 Bounding_box::Bounding_box(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){
     parameters[0] = xmin;
     parameters[1] = xmax;
     parameters[2] = ymin;
     parameters[3] = ymax;
     parameters[4] = zmin;
     parameters[5] = zmax;
 }

bool Bounding_box::ray_box_intersection(Ray R){
    float tmin_y, tmax_y, tmin_z, tmax_z;
    vector3 inv_direction(1/(R.direction.x()), 1/R.direction.y(), 1/R.direction.z());
    tmin = (parameters[0]- R.origin.x())*inv_direction.x();
    tmax = (parameters[1]-R.origin.x())*inv_direction.x();
    tmin_y = (parameters[2]- R.origin.y())*inv_direction.y();
    tmax_y = (parameters[3]- R.origin.y())*inv_direction.y();

     if (tmin > tmax)std::swap(tmin, tmax); 
     if (tmin_y > tmax_y)std::swap(tmin_y, tmax_y);  
     if ((tmin > tmax_y)||(tmin_y>tmax)){
         return 0;
     }
     if (tmin_y > tmin){
        tmin = tmin_y;
     }
     if (tmax_y < tmax){
         tmax = tmax_y;
     }        
     tmin_z = (parameters[4]- R.origin.z())*inv_direction.z();
     tmax_z = (parameters[5]- R.origin.z())*inv_direction.z();
   
    //  if (tmin_z > tmax_z)std::swap(tmin_z, tmax_z);

    //  if ((tmin > tmax_z)||(tmin_z>tmax)){
    //      return 0;
    //  }
     return 1;     
 }

void search_tree::traverse_tree(search_tree*root, Ray R, std::vector<int> *output){
    if(((root->left_node!=nullptr))||((root->right_node!=nullptr))){
        Bounding_box B_right(root->right_node->parameters[0],root->right_node->parameters[1], root->right_node->parameters[2],root->right_node->parameters[3],root->right_node->parameters[4],root->right_node->parameters[5]);
        Bounding_box B_left(root->left_node->parameters[0],root->left_node->parameters[1], root->left_node->parameters[2],root->left_node->parameters[3],root->left_node->parameters[4],root->left_node->parameters[5]);
        if ((root->right_node!=nullptr)&&(B_right.ray_box_intersection(R)==1)) {
            traverse_tree(root->right_node,R, output);
        }  
        if ((root->left_node!=nullptr)&&(B_left.ray_box_intersection(R)==1)){
            traverse_tree(root->left_node, R, output);
        }   
    } 
    if(((root->left_node==nullptr))&&((root->right_node==nullptr))){
        Bounding_box B_root(root->parameters[0],root->parameters[1], root->parameters[2],root->parameters[3],root->parameters[4],root->parameters[5]);
        if((B_root.ray_box_intersection(R)==1)){        
            for (int i = 0; i<root->number_of_node_faces; i++){
                (*output).push_back( root->faces_in_node[i]);
            }         
        }   
    }
}

void search_tree::find_parameters(int i, float* vertices,int*faces, std::vector<float> *parameters, std::vector<float> initial_parameters){
    float xmin = initial_parameters[0], ymin = initial_parameters[2], zmin = initial_parameters[4], xmax=initial_parameters[1], ymax=initial_parameters[3], zmax = initial_parameters[5];
    for(int j=0; j<3; j++){
        if (vertices[3*(faces[3*i+j]-1)]< xmin){
        xmin = vertices[3*(faces[3*i+j]-1)];
        }
        if (vertices[3*(faces[3*i+j]-1)+1]< ymin){
        ymin = vertices[3*(faces[3*i+j]-1)+1];
        }
        if (vertices[3*(faces[3*i+j]-1)+2]< zmin){
        zmin = vertices[3*(faces[3*i+j]-1)+2];
        }
        if (vertices[3*(faces[3*i+j]-1)]> xmax){
        xmax = vertices[3*(faces[3*i+j]-1)];
        }
        if (vertices[3*(faces[3*i+j]-1)+1]> ymax){
        ymax =vertices[3*(faces[3*i+j]-1)+1];
        }
        if (vertices[3*(faces[3*i+j]-1)+2]> zmax){
        zmax = vertices[3*(faces[3*i+j]-1)+2];
        }
    }   
    *parameters = {xmin, xmax, ymin, ymax, zmin, zmax};

}

void search_tree::leaf_nodes(float* vertices, int*faces, int number_of_faces, std::vector<search_tree*> *leaf_nodes){
    std::vector<float> parameters;
    std::vector<float> initial_parameters = {infinity, -1*infinity, infinity, -1*infinity, infinity, -1*infinity};
    for(int i = 0; i<number_of_faces; i++){ // make list of leaf_nodes.
        search_tree* leaf = new search_tree;
        leaf->number_of_node_faces = 1; 
        leaf->faces_in_node= new int [1];
        leaf->faces_in_node[0]=i;
        parameters.clear();
        search_tree::find_parameters(i, vertices, faces, &parameters, initial_parameters);
        for(int j =0; j<6;j++){
            leaf->parameters[j] = parameters[j];
        }
        leaf->left_node=nullptr;
        leaf->right_node=nullptr;
        (*leaf_nodes).push_back(leaf);
    }    
}

void search_tree::build_tree(float* vertices, int* faces, std::vector<search_tree*>* leaf_nodes, search_tree**root){
    std::vector<float> parameters;
    int it=0;
   
    while((*leaf_nodes).size()>1){
        it=it+1;
        parameters.clear();
        search_tree* temp = new search_tree;
        search_tree* first = (*leaf_nodes)[0];
        search_tree* second = (*leaf_nodes)[1];
        temp->left_node =(*leaf_nodes)[0];
        temp->right_node = (*leaf_nodes)[1];
        temp->number_of_node_faces = first->number_of_node_faces+second->number_of_node_faces;
        (*leaf_nodes).erase ((*leaf_nodes).begin(),(*leaf_nodes).begin()+2);
        temp->faces_in_node = new int [temp->number_of_node_faces];
        for(int i=0; i<first->number_of_node_faces; i++){
            temp->faces_in_node[i] = first->faces_in_node[i];
        }
        for(int i = first->number_of_node_faces; i<first->number_of_node_faces+second->number_of_node_faces;i++){
            temp->faces_in_node[i] = second->faces_in_node[i-first->number_of_node_faces];
        }
        std::vector<float> initial_parameters = {infinity, -1*infinity, infinity, -1*infinity, infinity, -1*infinity};
        for(int i =0; i<temp->number_of_node_faces; i++){
            search_tree::find_parameters(temp->faces_in_node[i], vertices, faces, &parameters, initial_parameters);
            initial_parameters = {parameters[0],parameters[1],parameters[2], parameters[3], parameters[4], parameters[5]};              
        }
        for(int j =0; j<6;j++){
             temp->parameters[j] = parameters[j];
        }
        (*leaf_nodes).push_back(temp);
    }
       *root = (*leaf_nodes)[0];
        return;
}

void search_tree::delete_tree(search_tree* root){
    if((root->left_node)!=nullptr){
        delete_tree((root)->left_node);
    }
     if((root->right_node)!=nullptr){
        delete_tree((root)->right_node);
     }
    delete root->faces_in_node;
    delete root;
    root = nullptr;   
}