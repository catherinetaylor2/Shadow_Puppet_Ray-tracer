//---------------- Shadow Puppet Monte-Carlo Ray-tracer -----------------------------
//
// Created by Catherine Taylor
//
// Began May 2017
//
//Produces hard shadow puppets

#include <cstdlib>
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

    unsigned char * TextureData, * PuppetTexture;
	int ScreenTextureWidth, ScreenTextureHeight,  PuppetTextureWidth, PuppetTextureHeight;
	TextureData = readBMP("Textures/sheet_5.bmp", &ScreenTextureWidth, &ScreenTextureHeight); //Screen texture input
    if(TextureData == 0){
        std::cerr<<"Error: Screen texture does not exist \n";
        return -1;
    }
	PuppetTexture = readBMP("Textures/turtle_texture.bmp", &PuppetTextureWidth, &PuppetTextureHeight); //pupet texture input
    if(PuppetTexture == 0){
        std::cerr<<"Error: Puppet texture does not exist \n";
        return -1;
    }
    std::cout<<"Texture bitmaps loaded\n";

    int width, height;
	if(argc>1){
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
	else{
		width = 500;
		height = 500;
	}

    //Quad mesh inputs
    float *vertices, *normals, *Textures;
    int numberOfFaces, *faceVertices, *faceNormals, *faceTextures;
    ObjFile DinoMesh("Objects/quad.obj");
    if(DinoMesh.doesExist()==false){
        std::cerr<<"Error: Object does not exist \n";
        return -1;
    }
	DinoMesh.get_mesh_data(DinoMesh, &faceVertices, &faceNormals, &faceTextures, &Textures, &normals, &vertices, &numberOfFaces);
    search_tree* root; 
    std::vector<search_tree*> LeafNodes;
    search_tree::leaf_nodes(vertices, faceVertices, numberOfFaces, &LeafNodes);
	search_tree::build_tree(vertices, faceVertices, &LeafNodes, &root);
	std::cout<<"tree built \n";

    //Set up camera position
    vector3 eye(0.0f,0.0f,-75.0f), lookat(0.0f,0.0f,1.0f), lookup(0.0f,1.0f,-30.0f);

    //Set up scene and light position.
    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup);
    float LightLength = 0.25f;
    vector3 LightCentre(0.0f, 0.0f, 50.0f);
    light myLight(LightLength, 1.0f, LightCentre);

    int iterations = 50; //number of rays per pixel
	unsigned char *img = new unsigned char[3*myScene.get_x_res()*myScene.get_y_res()];

    for (int x = 0; x<3*myScene.get_x_res()*myScene.get_y_res(); x+=3){ //loops over all pixels
      
        int i, j;
        i=(x/(3))%(myScene.get_x_res());
        j=(x/(3))/(myScene.get_x_res());

        vector3 pixelCoord = vector3::add3(myScene.get_corner(), vector3::ScalarMultiply(1*i*myScene.get_ratio(),myScene.u()), vector3::ScalarMultiply(-1*j*myScene.get_ratio(),myScene.v()) ); //pixel poPointOnLighttion in world space.

        float value = 0.0f, PixelColourSum= 0.0f;
        int adaptive = 0, testIterations = 25 ; //initial values for adaptive sampling
        float* intersectionColours = new float[testIterations];

        #pragma omp parallel for
        for(int z =0; z <testIterations; ++z){
            vector3 PointOnLight = myLight.PointOnSource();
            vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
            rayDirections.normalize();
            vector3 LightRayDirection = vector3::subtract(pixelCoord, PointOnLight);
            LightRayDirection.normalize();
            Ray RayFired(pixelCoord, rayDirections);
            
            value += triangle::intersection_value(RayFired, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, myLight.get_normal(), LightRayDirection, &intersectionColours, z );
            
        }

        for(int z = 0; z<testIterations; ++z){
            PixelColourSum+= intersectionColours[z]; 
        }
        for(int z = 0; z<testIterations; ++z){
            if(((intersectionColours)[z]/PixelColourSum>1.0f/(float)testIterations)&&(PixelColourSum>0)){ //if one ray differs significantly then test more.
                adaptive = 1;
            }
        }    

        if(adaptive==1){ //if needed used adaptive and repeat above.
            #pragma omp parallel for 
            for (int l=0; l<iterations; ++l){
                
                vector3 PointOnLight = myLight.PointOnSource();
                vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
                rayDirections.normalize();
                vector3 LightRayDirection = vector3::subtract(pixelCoord, PointOnLight);
                LightRayDirection.normalize();
                Ray RayFired(pixelCoord, rayDirections);

                value += triangle::intersection_value(RayFired, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight,  myLight.get_normal(), LightRayDirection, &intersectionColours, 0 );
                
            }
        }
        delete[] intersectionColours;
        
        //Using Monte Carlo, average values.
        float R = TextureData[j*ScreenTextureWidth*3 + 3*i]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        float G = TextureData[j*ScreenTextureWidth*3 + 3*i+1]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        float B = TextureData[j*ScreenTextureWidth*3 + 3*i+2]*value/(float)(iterations*(adaptive==1)+testIterations)*myLight.get_illumination();
        
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
    delete [] TextureData;
    delete [] PuppetTexture;

    return 0;
}