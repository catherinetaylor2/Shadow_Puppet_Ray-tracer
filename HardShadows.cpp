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
#include "ReadObj.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include "BITMAP.hpp"
#include "scene.hpp"
#include "triangle.hpp"
#include "binarySearchTree.hpp"

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
    float *vertices, *normals, *Textures;     //Quad mesh inputs
    int numberOfFaces, *faceVertices, *faceNormals, *faceTextures, numberOfVertices;
    ObjFile DinoMesh("Objects/quad.obj");
    if(DinoMesh.doesExist()==false){
        std::cerr<<"Error: Object does not exist \n";
        return -1;
    }
	DinoMesh.getMeshData(DinoMesh, &faceVertices, &faceNormals, &faceTextures, &Textures, &normals, &vertices, &numberOfFaces, &numberOfVertices);
    binarySearchTree* root; 
    std::vector<binarySearchTree*> LeafNodes;
    binarySearchTree::findLeafNodes(vertices, faceVertices, numberOfFaces, &LeafNodes);
	binarySearchTree::buildTree(vertices, faceVertices, &LeafNodes, &root);
	std::cout<<"Inputs loaded \n";

    vector3 eye(0.0f,0.0f,-75.0f), lookat(0.0f,0.0f,1.0f), lookup(0.0f,1.0f,-30.0f);  //Set up camera position

    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup); //Set up scene and light position.
    float LightLength = 0.25f;
    vector3 LightCentre(0.0f, 0.0f, 50.0f);
    light myLight(LightLength, 1.0f, LightCentre);

    int iterations = 50; //number of rays per pixel
	unsigned char *img = new unsigned char[3*myScene.getXRes()*myScene.getYRes()];

    for (int x = 0; x<3*myScene.getXRes()*myScene.getYRes(); x+=3){ //loops over all pixels
      
        int i, j;
        i=(x/(3))%(myScene.getXRes());
        j=(x/(3))/(myScene.getXRes());

        vector3 pixelCoord = vector3::add3(myScene.getCorner(), vector3::ScalarMultiply(1*i*myScene.getRatio(),myScene.u()), vector3::ScalarMultiply(-1*j*myScene.getRatio(),myScene.v()) ); //pixel poPointOnLighttion in world space.

        int testIterations = 25 ; //initial values for adaptive sampling
        float value = 0.0f, PixelColourSum= 0.0f, *intersectionColours = new float[testIterations];
        bool adaptive = false;

        #pragma omp parallel for
        for(int z =0; z <testIterations; ++z){
            vector3 PointOnLight = myLight.PointOnSource();
            vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
            rayDirections.normalize();
            vector3 LightRayDirection = vector3::subtract(pixelCoord, PointOnLight);
            LightRayDirection.normalize();
            Ray RayFired(pixelCoord, rayDirections);
            value += triangle::getColour(RayFired, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, myLight.get_normal(), LightRayDirection, &intersectionColours, z );
        }

        for(int z = 0; z<testIterations; ++z){
            PixelColourSum += intersectionColours[z]; 
        }
        for(int z = 0; z<testIterations; ++z){
            if(((intersectionColours)[z]/PixelColourSum>1.0f/(float)testIterations)&&(PixelColourSum>0)){ //if one ray differs significantly then test more.
                adaptive = true;
            }
        }    

        if(adaptive==true){ //if needed used adaptive and repeat above.
            #pragma omp parallel for 
            for (int l=0; l<iterations; ++l){
                vector3 PointOnLight = myLight.PointOnSource();
                vector3 rayDirections = vector3::subtract(PointOnLight, pixelCoord); //from screen to light source
                rayDirections.normalize();
                vector3 LightRayDirection = vector3::subtract(pixelCoord, PointOnLight);
                LightRayDirection.normalize();
                Ray RayFired(pixelCoord, rayDirections);
                value += triangle::getColour(RayFired, root, vertices, faceVertices, faceTextures, Textures, PuppetTexture, PuppetTextureWidth, PuppetTextureHeight,  myLight.get_normal(), LightRayDirection, &intersectionColours, 0 );
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
    BITMAP_File_Header fileHeader;
    BITMAP_Info_Header infoHeader;
    fillBMPHeaders(&fileHeader, &infoHeader,  width, height);
    writeBMP(&fileHeader, &infoHeader,&image);
    for(auto x = height-1; x>=0; x--){
        for (auto y = 0; y < width; y++) {
            for(auto z =2; z>=0; z--){
            image<<img[x*width*3 + y*3+ z];
            }
        }
    }
    image.close();

	ObjFile::cleanUp(vertices, normals, Textures, faceVertices, faceNormals, faceTextures); //Clear up files.
    binarySearchTree::deleteTree(root);
    delete [] img;
    delete [] TextureData;
    delete [] PuppetTexture;

    return 0;
}