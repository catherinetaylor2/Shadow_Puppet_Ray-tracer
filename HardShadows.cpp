//---------------- Shadow Puppet Monte-Carlo Ray-tracer -----------------------------
//
// Created by Catherine Taylor
//
// Began May 2017
//
//Produces hard shadow puppets

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

int main(int argc, char* argv[] ){

//Screen texture input
    unsigned char * data; 
	int ScreenTextureWidth, ScreenTextureHeight;
	data = readBMP("Textures/sheet_5.bmp", &ScreenTextureWidth, &ScreenTextureHeight);
    std::cout<<"width "<<ScreenTextureWidth<<" "<<ScreenTextureHeight<<"\n";

//Puppet texture input
    unsigned char * PuppetTexture;
	int PuppetTextureWidth, PuppetTextureHeight;
	PuppetTexture = readBMP("Textures/dino_texture.bmp", &PuppetTextureWidth, &PuppetTextureHeight);
    std::cout<<"quad texture width "<<PuppetTextureWidth<<"  "<<PuppetTextureHeight<<"\n";

    int width, height;
	if(argc>1){
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
	else{
		width=1920;
		height=1080;
	}


for (int ObjFileInput = 1;  ObjFileInput<2; ObjFileInput++){
    std::string j;

        if(ObjFileInput < 10){
			j = "000"+std::to_string(ObjFileInput);
		}
		else if ((ObjFileInput>=10)&&(ObjFileInput<100)){
			j= "00" + std::to_string(ObjFileInput);
		}
		else{
			j= "0" + std::to_string(ObjFileInput);
		}

//Quad mesh inputs
    float *vertices, *normals, *Textures;
    int numberOfFaces, *faceVertices, *faceNormals, *faceTextures;
    ObjFile mesh_dino("Objects/quad.obj");
	mesh_dino.get_mesh_data(mesh_dino, &faceVertices, &faceNormals, &faceTextures, &Textures, &normals, &vertices, &numberOfFaces);
    search_tree* root; 
    std::vector<search_tree*> LeafaceNormalsodes;
    search_tree::leaf_nodes(vertices, faceVertices, numberOfFaces, &LeafaceNormalsodes);
	search_tree::build_tree(vertices, faceVertices, &LeafaceNormalsodes, &root);
	std::cout<<"tree built \n";

//Set up camera poPointOnLighttion
    vector3 eye(0.0f,0.0f,-75.0f), lookat(0.0f,0.0f,1.0f), lookup(0.0f,1.0f,-30.0f);

//Set up scene and light poPointOnLighttion.
    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup);
    float LightLength = 0.25f;
    vector3 LightCentre(0.0f, 0.0f, 50.0f);
    light myLight(LightLength, 1.0f, LightCentre);

    int iterations=10; //number of rays per pixel
	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];

    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){ //loops over all pixels
        bool viPointOnLightbility;
        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 pixelCoord = vector3::add3(myScene.get_corner(), vector3::vec_scal_mult(1*i*myScene.get_ratio(),myScene.get_u()), vector3::vec_scal_mult(-1*j*myScene.get_ratio(),myScene.get_v()) ); //pixel poPointOnLighttion in world space.

        float value = 0.0f, sum = 0.0f;
        int adaptive = 0, testIterations = 25 ; //initial values for adaptive sampling
        float* colours = new float[testIterations];

        #pragma omp parallel for
        for(int z =0; z <testIterations; z++){
            vector3 PointOnLight = myLight.point_on_source();
            vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
            rayDirections.normalize();
            vector3 L = vector3::subtract(pixelCoord, PointOnLight);
            L.normalize();
            Ray R(pixelCoord, rayDirections);
            
            float temp_value = triangle::intersection_value(R, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, myLight.get_normal(), L, &colours, z );
            value = value + temp_value;
        }

        for(int z = 0; z<testIterations; z++){
            sum += colours[z]; 
        }
        for(int z = 0; z<testIterations; z++){
            if(((colours)[z]/sum >1.0f/(float)testIterations)&&(sum>0)){ //if one ray differs PointOnLightgnificantly then test more.
                adaptive = 1;
            }
        }    

        if(adaptive==1){ //if needed used adaptive and repeat above.
            #pragma omp parallel for 
            for (int l=0; l<iterations;l++){
                vector3 PointOnLight = myLight.point_on_source();
                vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
                rayDirections.normalize();
                vector3 L = vector3::subtract(pixelCoord, PointOnLight);
                L.normalize();
                Ray R(pixelCoord, rayDirections);

                float temp_value = triangle::intersection_value(R, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight,  myLight.get_normal(), L, &colours, 0 );
                value = value + temp_value;
            }
        }
        delete[] colours;
        
//UPointOnLightng Monte Carlo, average values.
        float R = data[j*ScreenTextureWidth*3 + 3*i]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        float G = data[j*ScreenTextureWidth*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        float B = data[j*ScreenTextureWidth*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        
        if(R>255.0f){ //clamp at 255
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
  
    std::ofstream image("puppet.bmp", std::ios::out| std::ios::binary); //write to bmp file.
    BITMAP_File_Header file_header;
    BITMAP_Info_Header info_header;
    fill_bitmap_headers(&file_header, &info_header,  width, height);
    write_bitmap (&file_header, &info_header,&image);
    for(auto x = height-1; x>=0; x--){
        for (auto y = 0; y < width; y++) {
            for(auto z =2; z>=0; z--){
            image<<img[x*width*3 + y*3+ z];
            }
        }
    }
    image.close();



//Clear up files.
	ObjFile::clean_up(vertices, normals, Textures, faceVertices, faceNormals, faceTextures);
    search_tree::delete_tree(root);
    delete [] img;
}
    delete [] data;
    delete [] PuppetTexture;

    return 0;
}