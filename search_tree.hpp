#ifndef search_tree_hpp
#define search_tree_hpp

#include <iostream>
#include "vec3.hpp"
#include <vector>

class search_tree{
    public:
        search_tree *left_node;
        search_tree *right_node;
        float parameters [6];
        int*faces_in_node, number_of_node_faces;
        static void build_tree(float* vertices, int* faces, std::vector<search_tree*> *leaf_nodes, search_tree**root);
        static void traverse_tree(search_tree*root, vector3 eye, vector3 d, std::vector<int> *output);
        static void leaf_nodes(float* vertices, int*faces, int number_of_faces, std::vector<search_tree*> *leaf_nodes);
        static void find_parameters(int i, float* vertices, int*faces,  std::vector<float> *parameters, std::vector<float> initial_parameters);
        static void delete_tree(search_tree* root);
        void get_tree_data(float* vertices, int* faces, std::vector<search_tree*> *leaf_nodes, search_tree**root);
    private:
};
class Bounding_box{
    public:
        Bounding_box(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);  
        int ray_box_intersection(vector3 ray_point, vector3 ray_direction);    
        float get_tmin(void){
            return tmin;
        }
        float get_tmax(void){
            return tmax;
        }
    private:
        float parameters [6];
        float tmin;
        float tmax;
};

#endif