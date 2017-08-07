//---------------- Shadow Puppet Monte-Carlo Ray-tracer -----------------------------
//
// Created by Catherine Taylor
//
// Began May 2017
//
//Produces sofaceTextures shadow puppets

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
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
    
    int width, height;
	if(argc>1){
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
	else{
		width = 500, height = 500;
	}

//Input texture data
    unsigned char * ScreenData, *PuppetTexture;
	int ScreenTextureWidth, ScreenTextureHeight, PuppetTextureWidth, PuppetTextureHeight;
	ScreenData = readBMP("Textures/sheet6.bmp", &ScreenTextureWidth, &ScreenTextureHeight);
     if(ScreenData == 0){
        std::cerr<<"Error: Screen texture does not exist \n";
        return -1;
    }
	PuppetTexture = readBMP("Textures/seahorse_texture.bmp", &PuppetTextureWidth, &PuppetTextureHeight);
    if(PuppetTexture == 0){
        std::cerr<<"Error: Puppet textures does not exist \n";
        return -1;
    }

  

//Puppet mesh inputs
    float *vertices, *normals, *textures;
    int NumberOfFaces, *faceVertices, *faceNormals, *faceTextures;
    ObjFile PuppetMesh("Objects/quad.obj");
      if(PuppetMesh.doesExist()==false){
        std::cerr<<"Error: Object does not exist \n";
        return -1;
    }
	PuppetMesh.get_mesh_data(PuppetMesh, &faceVertices, &faceNormals, &faceTextures, &textures, &normals, &vertices, &NumberOfFaces);
    search_tree* root; 
    std::vector<search_tree*> leafNodes;
    search_tree::leaf_nodes(vertices, faceVertices, NumberOfFaces, &leafNodes);
	search_tree::build_tree(vertices, faceVertices, &leafNodes, &root);
	std::cout<<"tree built \n";

//Camera input
    vector3 eye(0.0f,0.0f,-75.0f), lookat(0.0f,0.0f,1.0f), lookup(0.0f,1.0f,-30.0f);

//Scene set up
    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup);
    float innerLightLength = 0.1f, outerLightLength = 8.0f;
    vector3 innerCentre(0.0f, 0.0f, 50.0f), outerCentre(0.0f, 0.0f, 70.0f);
    light innerLight(innerLightLength,1.0f, innerCentre);
    light outerLight(outerLightLength, 0.9f, outerCentre);

    int iterations=100;

	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];
    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){

        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 s = vector3::add3(myScene.get_corner(), vector3::ScalarMultiply(1*i*myScene.get_ratio(),myScene.u()), vector3::ScalarMultiply(-1*j*myScene.get_ratio(),myScene.v()) );

        float value = 0,sum =0;
        int adaptive = 0, testIterations = 25 ; 
        float* colours = new float[testIterations];
        #pragma omp parallel for
        for(int z =0; z <testIterations; z++){
            vector3 Si = innerLight.PointOnSource();
            vector3 rayDirection(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
            vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
            L.normalize();
            rayDirection.normalize();
            Ray R(s, rayDirection);

            vector3 Sil =outerLight.PointOnSource();
            vector3 rayDirectionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
            vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
            Ll.normalize();
            rayDirectionl.normalize();
            Ray Rl(s, rayDirectionl);

            float tempValue = triangle::intersection_value_s(Rl, R, root, vertices, faceVertices, faceTextures, textures,PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, innerLight.get_normal(), L, &colours,  z, innerLight, outerLight);        
            value = value + tempValue;
        }

        for(int z = 0; z<testIterations; z++){
            sum += colours[z]; 
        }
        for(int z = 0; z<testIterations; z++){
            if(((colours)[z]/sum >0.05)&&(sum>0)){
                adaptive = 1;
            }
        }
        delete [] colours;      

        if(adaptive==1){
            #pragma omp parallel for 
            for (int l=0; l<iterations;l++){
                vector3 Si = innerLight.PointOnSource();
                vector3 rayDirection(Si.x()-s.x(), Si.y()-s.y(), Si.z()-s.z());
                vector3 L(s.x() - Si.x(), s.y() - Si.y(), s.z()-Si.z());
                L.normalize();
                rayDirection.normalize();
                Ray R(s, rayDirection);

                vector3 Sil =outerLight.PointOnSource();
                vector3 rayDirectionl(Sil.x()-s.x(), Sil.y()-s.y(), Sil.z()-s.z());
                vector3 Ll(s.x() - Sil.x(), s.y() - Sil.y(), s.z()-Sil.z());
                Ll.normalize();
                rayDirectionl.normalize();
                Ray Rl(s, rayDirectionl);

                float tempValue =  triangle::intersection_value_s(Rl, R, root, vertices, faceVertices, faceTextures, textures,PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, innerLight.get_normal(), L, &colours,  0, innerLight, outerLight);        
                value = value + tempValue;
            }
        }

         //Spherical light ScreenData:
        float R = ScreenData[j*ScreenTextureWidth*3 + 3*i]*value/(float)(iterations*(adaptive==1)+testIterations);
        float G = ScreenData[j*ScreenTextureWidth*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+testIterations);
        float B = ScreenData[j*ScreenTextureWidth*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+testIterations);

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

	ObjFile::clean_up(vertices, normals, textures, faceVertices, faceNormals, faceTextures);
    search_tree::delete_tree(root);
    delete[]PuppetTexture;
    delete[] ScreenData;
    delete[] img;

    return 0;
}