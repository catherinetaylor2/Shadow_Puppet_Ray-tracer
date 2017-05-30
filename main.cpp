//---------------- Shadow Puppet Monte-Carlo Ray-tracer -----------------------------
//
// Created by Catherine Taylor
//
// Began May 2017

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <thread>
#include "Read_Obj.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include "BITMAP.hpp"
#include "scene.hpp"
#include "triangle.hpp"
#include "search_tree.hpp"

#define infinity FLT_MAX;

double uniform_random_number(void){
    return rand()/double(RAND_MAX);
}

int main(int argc, char* argv[] ){

    int width, height;
	if(argc>1){
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
	else{
		width=1920;
		height=1080;
	}

//Puppet mesh inputs
    float *V, *N, *VT;
    int F, *FV, *FN, *FT;
    ObjFile mesh_dino("dino_puppet_handle.obj");
	mesh_dino.get_mesh_data(mesh_dino, &FV, &FN, &FT, &VT, &N, &V, &F);
    search_tree* root; 
    std::vector<search_tree*> leaf_nodes;
    search_tree::leaf_nodes(V, FV, F, &leaf_nodes);
	search_tree::build_tree(V, FV, &leaf_nodes, &root);
	std::cout<<"tree built \n";

//light mesh
    float *V_l, *N_l, *VT_l;
    int F_l, *FV_l, *FN_l, *FT_l;
    ObjFile mesh_light("light_source.obj");
	mesh_dino.get_mesh_data(mesh_light, &FV_l, &FN_l, &FT_l, &VT_l, &N_l, &V_l, &F_l);
    search_tree* root_l; 
    std::vector<search_tree*> leaf_nodes_l;
    search_tree::leaf_nodes(V_l, FV_l, F_l, &leaf_nodes_l);
	search_tree::build_tree(V_l, FV_l, &leaf_nodes_l, &root_l);
	std::cout<<"tree built \n";


    vector3 eye(0.0f,0.0f,-30.0f);
    vector3 lookat(0.0f,0.0f,1.0f);
    vector3 lookup(0.0f,1.0f,-30.0f);

    scene myScene(width, height, 90.0f, 40.0f, eye, lookat, lookup);
    float light_length = 10.0f;
    light myLight(-light_length, 75.0f, 1.0f);
    vector3 centre = myLight.get_centre();
    vector3 tangent_v(0,1,0);
    vector3 tangent_u(1,0,0);
    vector3 plane_n(0,0,1);
    int iterations=100;
    
	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];
    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){
        bool visibility;
        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 s = vector3::vec_add3(myScene.get_corner(), vector3::vec_scal_mult(1*i*myScene.get_ratio(),myScene.get_u()), vector3::vec_scal_mult(-1*j*myScene.get_ratio(),myScene.get_v()) );
        
        // vector3 ray_direction(s.x() - eye.x(), s.y()-eye.y(), s.z()-eye.z());
        // ray_direction.normalize();
        // Ray R(eye, ray_direction);
        // int min_value = -1, *k ; 
        // float t_min = triangle::intersection_point(root_l, V_l, R,FV_l, &min_value, &k); 
        // if(min_value!=-1){
        //     I=1;            
        // }
        // else{
        //     I=0;
        //   //  std::cout<<"oi oi \n";
        // }
        // delete k;
        


        float value = 0, alpha = 0.2f;
        #pragma omp parallel for 
        for (int l=0; l<iterations;l++){
            float a = uniform_random_number();
            float b = uniform_random_number();
            vector3 Si = vector3::vec_add3(centre, vector3::vec_scal_mult((0.5 - a)*light_length,tangent_u), vector3::vec_scal_mult((0.5 -b)*light_length,tangent_v));
            vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
            
            ray_direction.normalize();
            Ray R(s, ray_direction);

            int min_value = -1, *k ;
            float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k);
            if(min_value !=-1){
                #pragma omp critical
                value= value+alpha;
            }
            else{
                #pragma omp critical
                value = value+1;
            }
            
            delete[] k;

        }
        vector3 ray_direction(centre.x() - s.x(), centre.y()-s.y(), centre.z()-s.z());
        ray_direction.normalize();

    //  if(I==0){
        img[x]= value*190.0f/iterations;
        img[x+1]= value*120.0f/iterations;
        img[x+2]= value*45.0f/iterations;
        // }
        // else{
        //     img[x]= value*255.0f/iterations;
        //     img[x+1]= I*value*255.0f/iterations;
        //     img[x+2]= I*value*255.0f/iterations;
        // }

    }
    std::ofstream image2("puppet.bmp", std::ios::out| std::ios::binary); 
    BITMAP_File_Header file_header;
    BITMAP_Info_Header info_header;
    fill_bitmap_headers(&file_header, &info_header,  width, height);
    write_bitmap (&file_header, &info_header,&image2);
    for(auto x = height-1; x>=0; x--){
        for (auto y = 0; y < width; y++) {
            for(auto z =2; z>=0; z--){
            image2<<img[x*width*3 + y*3+ z];
            }
        }
    }
    image2.close();

	ObjFile::clean_up(V,N, VT, FV, FN, FT);
    search_tree::delete_tree(root);
    ObjFile::clean_up(V_l,N_l, VT_l, FV_l, FN_l, FT_l);
    search_tree::delete_tree(root_l);
    delete [] img;

    return 0;
}