//---------------- Shadow Puppet Monte-Carlo Ray-tracer -----------------------------
//
// Created by Catherine Taylor
//
// Began May 2017
//
//Produces soft shadow puppets

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
#define PI 3.141592654f

int main(int argc, char* argv[] ){

    unsigned char * data;
	int texture_width, texture_height;
	data = readBMP("sheet6.bmp", &texture_width, &texture_height);
    std::cout<<"width "<<texture_width<<" "<<texture_height<<"\n";

    unsigned char * dino_tex;
	int dino_width, dino_height;
	dino_tex = readBMP("dino_texture2.bmp", &dino_width, &dino_height);

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
    ObjFile mesh_dino("quad2.obj");
	mesh_dino.get_mesh_data(mesh_dino, &FV, &FN, &FT, &VT, &N, &V, &F);
    search_tree* root; 
    std::vector<search_tree*> leaf_nodes;
    search_tree::leaf_nodes(V, FV, F, &leaf_nodes);
	search_tree::build_tree(V, FV, &leaf_nodes, &root);
	std::cout<<"tree built \n";

    int* edges;
    mesh_dino.get_boundary_edges(FV, &edges, F);

    vector3 eye(0.0f,0.0f,-75.0f);
    vector3 lookat(0.0f,0.0f,1.0f);
    vector3 lookup(0.0f,1.0f,-30.0f);

    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup);
    float inner_light_length = 0.1f,I;
    light inner_light(inner_light_length, 50.0f, 1.0f);
    float outer_light_length = 7.5f;
    light outer_light(outer_light_length, 70.0f,1.0f);

    vector3 c(0.0f, 45.0f,50.0f);
    sphere_light sphereLight(c, 4.0f);

    vector3 plane_n(0,0,-1);
    vector3 puppet(0,0,1);
    int iterations=100;

	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];
    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){
        bool visibility;
        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 s = vector3::vec_add3(myScene.get_corner(), vector3::vec_scal_mult(1*i*myScene.get_ratio(),myScene.get_u()), vector3::vec_scal_mult(-1*j*myScene.get_ratio(),myScene.get_v()) );

        float value = 0,sum =0;
        float value_rgb = 0;
        int adaptive = 0, test_iterations = 25 ; 
        float* colours = new float[test_iterations];
        #pragma omp parallel for
        for(int z =0; z <test_iterations; z++){
            vector3 Si = inner_light.point_on_source();
            vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
            vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
            L.normalize();
            ray_direction.normalize();
            Ray R(s, ray_direction);

            vector3 Sil =outer_light.point_on_source();
            vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
            vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            Ll.normalize();
            ray_directionl.normalize();
            Ray Rl(s, ray_directionl);

            int min_value = -1, *k;
            float t_min = triangle::intersection_point(root, V, Rl,FV, &min_value, &k);
           
            int min_value2 = -1, *k2;
            float t_min2 = triangle::intersection_point(root, V, R,FV, &min_value2, &k2);            

            if(min_value2!=-1){              
                if(min_value!=-1){
                    int triangle = k[min_value+1];
                    float* colour = new float[3];
                    triangle::get_texture_value(triangle, FV, V, Rl, dino_tex, FT, VT, dino_width, dino_height, &colour); 

                    if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){
                        #pragma omp critical
                        value = value+0.2f;
                    }
                    else{
                        #pragma omp critical
                        value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                    }  
                    delete[] colour;
                    delete[] k;
                }
                else{
                    int triangle = k2[min_value2+1];
                    float* colour = new float[3];
                    triangle::get_texture_value(triangle, FV, V, R, dino_tex, FT, VT, dino_width, dino_height, &colour); 
                    if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){
                        #pragma omp critical
                        value = value+1.0f;                        
                    }
                    else{
                        #pragma omp critical
                        value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                    }  
                    delete[] colour;
                    delete[] k2;
                }
            }
            else{
                #pragma omp critical
                value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                #pragma omp critical
                value_rgb = value_rgb ;              
            }
            colours[z]=value; 
        }

        for(int z = 0; z<test_iterations; z++){
            sum += colours[z]; 
        }
        for(int z = 0; z<test_iterations; z++){
            if(((colours)[z]/sum >0.05)&&(sum>0)){
                adaptive = 1;
            }
        }
        delete [] colours;      

        if(adaptive==1){
            #pragma omp parallel for 
            for (int l=0; l<iterations;l++){
                vector3 Si = inner_light.point_on_source();
                vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
                vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
                L.normalize();
                ray_direction.normalize();
                Ray R(s, ray_direction);

                vector3 Sil =outer_light.point_on_source();
                vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
                vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
                Ll.normalize();
                ray_directionl.normalize();
                Ray Rl(s, ray_directionl);

                int min_value = -1, *k , min_value2 = -1, *k2 ;
                float t_min = triangle::intersection_point(root, V, Rl,FV, &min_value, &k);
                float t_min2 = triangle::intersection_point(root, V, R,FV, &min_value2, &k2);       
                
                if(min_value2!=-1){              
                    if(min_value!=-1){
                        int triangle = k[min_value+1];
                        float* colour = new float[3];
                        triangle::get_texture_value(triangle, FV, V, Rl, dino_tex, FT, VT, dino_width, dino_height, &colour); 

                        if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){
                            #pragma omp critical
                            value = value+0.2f;
                        }
                        else{
                            #pragma omp critical
                            value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                        }  
                        delete[] colour;
                        delete[] k;
                    }
                    else{
                        int triangle = k2[min_value2+1];
                        float* colour = new float[3];
                        triangle::get_texture_value(triangle, FV, V, R, dino_tex, FT, VT, dino_width, dino_height, &colour); 
                        if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){
                            #pragma omp critical
                            value = value+1.0f;     
                        }
                        else{
                            #pragma omp critical
                            value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                        }  
                        delete[] colour;
                        delete[] k2;
                    }
                }
                else{
                    #pragma omp critical
                    value = value+1.3f*pow(vector3::dotproduct(plane_n, L),10.0f);
                    #pragma omp critical
                    value_rgb = value_rgb ;              
                }
            }
         }

         //Spherical light data:
        float R = data[j*texture_width*3 + 3*i]*value/(float)(iterations*(adaptive==1)+test_iterations);
        float G = data[j*texture_width*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+test_iterations);//*221.0f/255.0f;
        float B = data[j*texture_width*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+test_iterations);//*204.0f/255.0f;

        if(R>255.0f){
            R = 255.0f;
        }
        if(G>255.0f){
            G=255.0f;
        }
        if(B>255.0f){
            B=255.0f;
        }
        img[x]=R;
        img[x+1]= G;
        img[x+2]=B ;  
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

    std::cout<<"done \n";

	ObjFile::clean_up(V,N, VT, FV, FN, FT);
    search_tree::delete_tree(root);
    delete[] dino_tex;
    delete[] data;
    delete[] edges;
    delete[] img;

    return 0;
}