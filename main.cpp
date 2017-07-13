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
#include <algorithm>
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
	data = readBMP("sheet_5.bmp", &texture_width, &texture_height);
    std::cout<<"width "<<texture_width<<" "<<texture_height<<"\n";

    unsigned char * dino_tex;
	int dino_width, dino_height;
	dino_tex = readBMP("dino_texture2.bmp", &dino_width, &dino_height);
    std::cout<<"w "<<dino_width<<" "<<dino_height<<"\n";

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
    ObjFile mesh_dino("quad.obj");
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
    float light_length = 0.25f,I;
    light myLight(light_length, 50.0f, 1.0f);

    vector3 plane_n(0,0,-1);
    vector3 puppet(0,0,1);
    

    vector3 V1(myLight.get_xmin(), myLight.get_ymin(), myLight.get_z());
    vector3 V2(myLight.get_xmin(), myLight.get_ymax(), myLight.get_z());
    vector3 V3(myLight.get_xmax(), myLight.get_ymin(), myLight.get_z());
    vector3 V4(myLight.get_xmax(), myLight.get_ymax(), myLight.get_z());
    triangle light_upper(V1, V4, V2);
    triangle light_lower(V3, V4, V1);

    int iterations=1;
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


            int min_value = -1, *k;
            float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k);
            
                    
            if(min_value!=-1){

                int triangle = k[min_value+1];
                float* colour = new float[3];
               
               triangle::get_texture_value(triangle, FV, V, R, dino_tex, FT, VT, dino_width, dino_height, &colour);      

   
            if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){

                    #pragma omp critical
                    value = value+0.2f;
                    }
                      else{
                    #pragma omp critical
                    value = value+1.3*pow(vector3::dotproduct(plane_n, L),50.0f)+0.4f;
                    #pragma omp critical
                    value_rgb = value_rgb+0.0f;
                }  
           delete[] colour;
            }
         
           else{
               #pragma omp critical
                value = value+1.3*pow(vector3::dotproduct(plane_n, L),50.0f)+0.4f;
                #pragma omp critical
                value_rgb = value_rgb ;              
           }
           colours[z]=value; 
           delete[] k;
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

//         if(adaptive==1){
//        std::cout<<"line 224 \n";
//          //  #pragma omp parallel for 
//             for (int l=0; l<iterations;l++){
//                 vector3 Si = myLight.point_on_source();
//                 vector3 ray_direction(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
//                 vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
//                 L.normalize();
//                 ray_direction.normalize();
//                 Ray R(s, ray_direction);

//                 int min_value = -1, *k ;
//                 float t_min = triangle::intersection_point(root, V, R,FV, &min_value, &k);
                
          
//                 if(min_value!=-1){
//                        vector3 POI = vector3::vec_add(s, vector3::vec_scal_mult(t_min, ray_direction));
//                                 // std::cout<<"s "<<s.x()<<" "<<s.y()<<" "<<s.z()<<"\n";
//                                 // std::cout<<"ray"<<ray_direction.x()<<" "<<ray_direction.y()<<" "<<ray_direction.z()<<"\n";
//                                 // std::cout<<"t "<<t_min<<"\n";

//    float denominator;
//     int triangle = k[min_value+1];
    
//    // std::cout<<"POI "<<POI.x()<<" "<<POI.y()<<" "<<POI.z()<<"\n";
//     int c_m1, c_m2, c_m3, n1, n2, n3;
//     c_m1 = FV[3*triangle] -1, c_m2 = FV[3*triangle+1]-1, c_m3 = FV[3*triangle+2] -1 ;
//    // n1 = FN[3*triangle]-1, n2 = FN[3*triangle+1]-1, n3= FN[3*triangle+2]-1;
//     vector3 point1(V[3*c_m1], V[3*c_m1+1], V[3*c_m1+2]);
//     vector3 point2(V[3*c_m2], V[3*c_m2+1], V[3*c_m2+2]);
//     vector3 point3(V[3*c_m3], V[3*c_m3+1], V[3*c_m3+2]);
//     // vector3 N1(N[3*n1], N[3*n1+1], N[3*n1+2]);
//     // vector3 N2(N[3*n2], N[3*n2+1], N[3*n2+2]);
//     // vector3 N3(N[3*n3], N[3*n3+1], N[3*n3+2]);
//     vector3 T = vector3::vec_add(R.get_origin(), vector3::vec_scal_mult(-1, point1));
//     vector3 E1 = vector3::vec_add(point2, vector3::vec_scal_mult(-1, point1));
//     vector3 E2 = vector3::vec_add(point3, vector3::vec_scal_mult(-1, point1));
//     denominator = vector3::dotproduct( vector3::crossproduct(R.get_direction(), E2), E1);

//     vector3 M (vector3::dotproduct(vector3::crossproduct(T, E1), E2),vector3::dotproduct(vector3::crossproduct(R.get_direction(), E2), T),vector3::dotproduct( vector3::crossproduct(T, E1), R.get_direction()));
//     vector3 tuv = vector3::vec_scal_mult(1.0f/(float)denominator, M);
//     // vector3 N = vector3::vec_add3(vector3::vec_scal_mult(std::max((float)(1-(tuv.y()+tuv.z())),0.0f),N1),vector3::vec_scal_mult(std::max(tuv.y(), 0.0f),N2),vector3::vec_scal_mult(std::max(tuv.z(),0.0f),N3));
//     float *barycentric = new float[3];
//     (barycentric)[0] = 1.0f-(tuv.y()+tuv.z());
//     (barycentric)[1] = tuv.y();
//     (barycentric)[2] = tuv.z();

//     int vt1 = FT[3*triangle]-1, vt2 = FT[3*triangle+1]-1, vt3 = FT[3*triangle+2]-1;
//     float vt_1x = VT[2*vt1], vt_1y = VT[2*vt1+1], vt_2x = VT[2*vt2], vt_2y = VT[2*vt2+1], vt_3x = VT[2*vt3], vt_3y = VT[2*vt3+1];

//      float u_coord, v_coord, alpha, beta, v12r, v12g, v12b, v34r, v34g, v34b;
//     int v1x,v1y, v2x, v4y;
//     u_coord = (barycentric[0]*vt_1x +barycentric[1]*vt_2x+barycentric[2]*vt_3x)*dino_width;
//     v_coord = (barycentric[0]*vt_1y +barycentric[1]*vt_2y+barycentric[2]*vt_3y)*dino_height;
//    // std::cout<<"line 146 "<<u_coord<<" "<<v_coord<<"\n";
//     v1x = (int)floor(u_coord);
//     v1y = (int)ceil(v_coord);
//     v2x = (int)ceil(u_coord);
//     v4y = (int)floor(v_coord);
//     if (v1x<0){
//         v1x=0;
//     }
//     if (v2x<0){
//         v2x=0;
//     }
//     if (v1y<0){
//         v1y=0;
//     }
//     if (v4y<0){
//         v4y=0;
//     }

//     alpha = (float)(u_coord - (v2x - v1x)*v1x)/(float) (v2x - v1x);
//     beta = (float)(v_coord - (v1y - v4y)*v4y)/(float) (v1y - v4y);

//     if (alpha >1){
//         alpha=1;
//     }
//     if(beta>1){
//         beta =1;
//     }
//  //std::cout<<"line 172 "<<v1y<<" "<<v1x<<"\n";
//     v12r = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x] +  alpha*dino_tex[v1y*dino_width*3 + 3*v2x];
// //std::cout<<"line 174 \n";
//     v12g = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x+1] +alpha*dino_tex[v1y*dino_width*3 + 3*v2x+1];
//     v12b = (1-alpha)*dino_tex[v1y*dino_width*3 + 3*v1x+2] +alpha*dino_tex[v1y*dino_width*3 + 3*v2x+2];
// //std::cout<<"line 176 \n";
//     v34g =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x+1] +  alpha*dino_tex[v4y*dino_width*3 + 3*v2x+1];
//     v34r =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x] +    alpha*dino_tex[v4y*dino_width*3 + 3*v2x];
//     v34b =  (1-alpha)*dino_tex[v4y*dino_width*3 + 3*v1x+2] +  alpha*dino_tex[v4y*dino_width*3 + 3*v2x+2];

//     float colour[] = {(1-beta)*v12r + beta*v34r, (1-beta)*v12g + beta*v34g, (1-beta)*v12b + beta*v34b} ;
// //std::cout<<"line 182 \n";
//                     if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){

//                     #pragma omp critical
//                     value = value+0.2f;
//                     }
//                       else{
//                     #pragma omp critical
//                     value = value+1.3*pow(vector3::dotproduct(plane_n, L),50.0f)+0.4f;
//                     #pragma omp critical
//                     value_rgb = value_rgb+0.0f;
//                 }  
//                 }
//                 else{
//                     #pragma omp critical
//                     value = value+1.3*pow(vector3::dotproduct(plane_n, L),50.0f)+0.4f;
//                     #pragma omp critical
//                     value_rgb = value_rgb+0.0f;
//                 }  
//                  delete[] k;    
//             }
              
//         }

        float R = data[j*texture_width*3 + 3*i]*value/(float)(iterations*(adaptive==1)+test_iterations);
        float G = data[j*texture_width*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+test_iterations);
        float B = data[j*texture_width*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+test_iterations);
  //  std::cout<<"line 342 \n";
        // float R = data[j*texture_width*3 + 3*i]*value/(float)(iterations*(adaptive==1)+test_iterations);//+0*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
        // float G = 221.0f/255.0f*data[j*texture_width*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+test_iterations);//+149*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
        // float B = 204.0f/255.0f*data[j*texture_width*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+test_iterations);//+30*value_rgb/(float)(iterations*(adaptive==1)+test_iterations);
        
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

std::cout<<"line 376 \n";
   delete [] img;
  
	ObjFile::clean_up(V,N, VT, FV, FN, FT);
    search_tree::delete_tree(root);
  //  delete [] img;
    delete [] data;
    delete [] edges;
    delete [] dino_tex;



    return 0;
}