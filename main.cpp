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
    ObjFile mesh_dino("dino_puppet_simple.obj");
	mesh_dino.get_mesh_data(mesh_dino, &FV, &FN, &FT, &VT, &N, &V, &F);
    search_tree* root; 
    std::vector<search_tree*> leaf_nodes;
    search_tree::leaf_nodes(V, FV, F, &leaf_nodes);
	search_tree::build_tree(V, FV, &leaf_nodes, &root);
	std::cout<<"tree built \n";


    vector3 eye(0.0f,0.0f,-30.0f);
    vector3 lookat(0.0f,0.0f,1.0f);
    vector3 lookup(0.0f,1.0f,-50.0f);

    scene myScene(width, height, 90.0f, 40.0f, eye, lookat, lookup);
    light myLight(-10.0f, 10.0f, -10.0f, 10.0f, 75.0f, 1.0f);
    vector3 centre = myLight.get_centre();
    vector3 tangent_v(0,1,0);
    vector3 tangent_u(1,0,0);
    float light_length = 10.0f, I =0;
    int iterations=10;
    
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
        // float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k); 
        // if(min_value !=-1){ 
        //     I=1;
        // }


        float value = 0;
        for (int l=0; l<iterations;l++){
            vector3 Si = vector3::vec_add3(centre, vector3::vec_scal_mult((0.5 - uniform_random_number())*light_length,tangent_u), vector3::vec_scal_mult((0.5 - uniform_random_number())*light_length,tangent_v));
        
            vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
            ray_direction.normalize();
            Ray R(s, ray_direction);

            int min_value = -1, *k ;
            float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k);
            if(min_value !=-1){
                value= value+0;
            }
            else{
                value = value+1;
            }
            
            delete k;

        }
       // if(I==0){
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
    delete [] img;

    return 0;
}