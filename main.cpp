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
    ObjFile mesh_dino("dino_puppet.obj");
	mesh_dino.get_mesh_data(mesh_dino, &FV, &FN, &FT, &VT, &N, &V, &F);
    search_tree* root; 
    std::vector<search_tree*> leaf_nodes;
    search_tree::leaf_nodes(V, FV, F, &leaf_nodes);
	search_tree::build_tree(V, FV, &leaf_nodes, &root);
	std::cout<<"tree built \n";

    vector3 eye(0.0f,0.0f,-50.0f);
    vector3 lookat(0.0f,0.0f,1.0f);
    vector3 lookup(0.0f,1.0f,-50.0f);
    vector3 light_centre(0.0f,0.0f,-150.0f);

    scene myScene(width, height, 90.0f, 3.0f, eye, lookat, lookup);
    light myLight(light_centre, 5.0f, 1.0f);
    
	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];
    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){
        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 s = vector3::vec_add3(myScene.get_corner(), vector3::vec_scal_mult(1*i*myScene.get_ratio(),myScene.get_u()), vector3::vec_scal_mult(-1*j*myScene.get_ratio(),myScene.get_v()) );
        vector3 ray_direction(s.x() - eye.x(), s.y()-eye.y(), s.z()-eye.z());
        Ray R(eye, ray_direction);

        // float t_min = 0;
        // for(int k =0; k<F; k++){
        //     int index1 = FV[3*k]-1, index2 = FV[3*k+1]-1, index3 = FV[3*k+2]-1;
        //     vector3 V1(V[3*index1], V[3*index1+1], V[3*index1+2]);
        //     vector3 V2(V[3*index2], V[3*index2+1], V[3*index2+2]);
        //     vector3 V3(V[3*index3], V[3*index3+1], V[3*index3+2]);
        //     triangle T(V1,V2,V3);
        //     bool t = T.ray_triangle_intersection(R);
        //     if(t!=0){
        //         t_min=1;
        //     }
        // }
        //     if(t_min!=0){
        //         img[x]=255;
        //     }
        //     else{
        //         img[x]=0;
        //     }
        //     img[x+1]=0;
        //     img[x+2]=0;

            
        int min_value = -1, *k ;
        float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k);
         if(min_value !=-1){
             int m = k[min_value+1];
            int index1 = FV[3*m]-1, index2 = FV[3*m+1]-1, index3 = FV[3*m+2]-1;
            vector3 V1(V[3*index1], V[3*index1+1], V[3*index1+2]);
            vector3 V2(V[3*index2], V[3*index2+1], V[3*index2+2]);
            vector3 V3(V[3*index3], V[3*index3+1], V[3*index3+2]);

            triangle T(V1,V2,V3);
            bool t = T.ray_triangle_intersection(R);
            if(t!=0){
                img[x]=255;
            }
            else{
                img[x]=0;
            }
        }
            else{
                img[x]=0;
            }
            img[x+1]=0;
            img[x+2]=0;
            delete k;
            
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