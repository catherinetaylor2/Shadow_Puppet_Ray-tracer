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
#define PI 3.141592654f

// double uniform_random_number(void){
//     return rand()/double(RAND_MAX);
// }

int main(int argc, char* argv[] ){

    unsigned char * data;
	int texture_width, texture_height;
	data = readBMP("sheet.bmp", &texture_width, &texture_height);
    std::cout<<"width "<<texture_width<<" "<<texture_height<<"\n";

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
    ObjFile mesh_dino("dino_puppet_handle2.obj");
	mesh_dino.get_mesh_data(mesh_dino, &FV, &FN, &FT, &VT, &N, &V, &F);
    search_tree* root; 
    std::vector<search_tree*> leaf_nodes;
    search_tree::leaf_nodes(V, FV, F, &leaf_nodes);
	search_tree::build_tree(V, FV, &leaf_nodes, &root);
	std::cout<<"tree built \n";

    vector3 eye(0.0f,0.0f,-30.0f);
    vector3 lookat(0.0f,0.0f,1.0f);
    vector3 lookup(0.0f,1.0f,-30.0f);

    scene myScene(width, height, 90.0f, 40.0f, eye, lookat, lookup);
    float light_length = 0.25f,I;
    light myLight(light_length, 100.0f, 1.0f);
    light light2(5.0f, 50.0f,1.0f);

    vector3 c(0.0f, 5.0f,50.0f);
    sphere_light sphereLight(c, 4.0f);

    vector3 plane_n(0,0,-1);
    int iterations=10;

    vector3 V1(myLight.get_xmin(), myLight.get_ymin(), myLight.get_z());
    vector3 V2(myLight.get_xmin(), myLight.get_ymax(), myLight.get_z());
    vector3 V3(myLight.get_xmax(), myLight.get_ymin(), myLight.get_z());
    vector3 V4(myLight.get_xmax(), myLight.get_ymax(), myLight.get_z());
    triangle light_upper(V1, V4, V2);
    triangle light_lower(V3, V4, V1);
     
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
            vector3 Si = myLight.point_on_source();
            vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
            vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
            L.normalize();
            ray_direction.normalize();
            Ray R(s, ray_direction);

            // vector3 Sil = sphereLight.point_on_source();
            // vector3 light_normal(Si.x()-sphereLight.get_centre().x(), Si.y()-sphereLight.get_centre().y(),Si.z()-sphereLight.get_centre().z());
            // light_normal.normalize();
            // vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            // Ll.normalize();
            // vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.z(), Sil.z()-s.z());
            // ray_directionl.normalize();
            // Ray Rl(s, ray_directionl);

            vector3 Sil =light2.point_on_source();
            vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
            vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            Ll.normalize();
            ray_directionl.normalize();
            Ray Rl(s, ray_directionl);


            int min_value = -1, *k;
            float t_min = triangle::intersection_point(root, V, Rl,FV, &min_value, &k);
            delete[] k;
            int min_value2 = -1, *k2;
            float t_min2 = triangle::intersection_point(root, V, R,FV, &min_value2, &k2);
            delete[] k2;

          if(min_value2!=-1){              
            if(min_value!=-1){
                    #pragma omp critical
                    value = value+0.2f;
            }
            else{
                value = value+1.2*pow(vector3::dotproduct(plane_n, L),30.0f);
            }
           }
           else{
               #pragma omp critical
                value = value+1.5*pow(vector3::dotproduct(plane_n, L),50.0f);
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
                vector3 Si = myLight.point_on_source();
                vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
                vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
                L.normalize();
                ray_direction.normalize();
                Ray R(s, ray_direction);


            // vector3 Sil = sphereLight.point_on_source();
            // // vector3 light_normal(Si.x()-sphereLight.get_centre().x(), Si.y()-sphereLight.get_centre().y(),Si.z()-sphereLight.get_centre().z());
            // // light_normal.normalize();
            // vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            // Ll.normalize();
            // vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.z(), Sil.z()-s.z());
            // ray_directionl.normalize();
            // Ray Rl(s, ray_directionl);

            vector3 Sil =light2.point_on_source();
            vector3 ray_directionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
            vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            Ll.normalize();
            ray_directionl.normalize();
            Ray Rl(s, ray_directionl);

            // vector3 ray_direction(c.x() - s.x(), c.y()-s.y(), c.z()-s.z());
            // ray_direction.normalize();
            // Ray R(s, ray_direction);

                int min_value = -1, *k , min_value2 = -1, *k2 ;
                float t_min = triangle::intersection_point(root, V, Rl,FV, &min_value, &k);
                delete[] k;
                float t_min2 = triangle::intersection_point(root, V, R,FV, &min_value2, &k2);
                delete[] k2;
        
                
                if(min_value2!=-1){
                    if(min_value!=-1){
                            #pragma omp critical
                            value = value+0.2f;
                        }
                        else{
                            value = value+1.2*pow(vector3::dotproduct(plane_n, L),30.0f);
                        }
                }
                else{
                    #pragma omp critical
                    value = value+1.5*pow(vector3::dotproduct(plane_n, L),50.0f) ;
                    #pragma omp critical
                    value_rgb = value_rgb+0.0f;
                }    
               
            }
        }


        float R = data[j*texture_width*3 + 3*i]*value/(float)(iterations*(adaptive==1)+test_iterations);//+0*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
        float G = data[j*texture_width*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+test_iterations);//+149*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
        float B = data[j*texture_width*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+test_iterations);//+30*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
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

	ObjFile::clean_up(V,N, VT, FV, FN, FT);
    search_tree::delete_tree(root);
    delete [] img;
    delete [] data;



    return 0;
}