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
#include <algorithm>
#include "Read_Obj.hpp"
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

    unsigned char * ScreenData, *PuppetTexture; //Input texture data
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

    float *vertices, *normals, *textures; //Puppet mesh inputs
    int NumberOfFaces, *faceVertices, *faceNormals, *faceTextures;
    ObjFile PuppetMesh("Objects/quad.obj");
    if(PuppetMesh.doesExist()==false){
        std::cerr<<"Error: Object does not exist \n";
        return -1;
    }
	PuppetMesh.getMeshData(PuppetMesh, &faceVertices, &faceNormals, &faceTextures, &textures, &normals, &vertices, &NumberOfFaces);
    binarySearchTree* root; 
    std::vector<binarySearchTree*> leafNodes;
    binarySearchTree::findLeafNodes(vertices, faceVertices, NumberOfFaces, &leafNodes);
	binarySearchTree::buildTree(vertices, faceVertices, &leafNodes, &root);
	std::cout<<"tree built \n";

    vector3 eye(0.0f,0.0f,-75.0f), lookat(0.0f,0.0f,1.0f), lookup(0.0f,1.0f,-30.0f); //Camera input

    scene myScene(width, height, 90.0f, 60.0f, eye, lookat, lookup); //Scene set up
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

        vector3 PointOnScreen = vector3::add3(myScene.get_corner(), vector3::ScalarMultiply(1*i*myScene.get_ratio(),myScene.u()), vector3::ScalarMultiply(-1*j*myScene.get_ratio(),myScene.v()) );

        int testIterations = 25 ; 
        float value = 0.0f, PixelColourSum = 0.0f, *intersectionColours = new float[testIterations];
        bool adaptive = false;

        #pragma omp parallel for
        for(int z =0; z <testIterations; z++){
            vector3 PointOnInnerLight = innerLight.PointOnSource();
            vector3 rayDirection(PointOnInnerLight.x()-PointOnScreen.x(), PointOnInnerLight.y()-PointOnScreen.y(), PointOnInnerLight.z()-PointOnScreen.z());
            vector3 LightDirInner(PointOnScreen.x() - PointOnInnerLight.x(), PointOnScreen.y() - PointOnInnerLight.y(), PointOnScreen.z()-PointOnInnerLight.z());
            LightDirInner.normalize();
            rayDirection.normalize();
            Ray rayInner(PointOnScreen, rayDirection);

            vector3 PointOnOuterLight = outerLight.PointOnSource();
            vector3 rayDirectionOuter(PointOnOuterLight.x()-PointOnScreen.x(), PointOnOuterLight.y()-PointOnScreen.y(), PointOnOuterLight.z()-PointOnScreen.z());
            vector3 LightDirOuter(PointOnScreen.x() - PointOnOuterLight.x(), PointOnScreen.y() - PointOnOuterLight.y(), PointOnScreen.z()-PointOnOuterLight.z());
            LightDirOuter.normalize();
            rayDirectionOuter.normalize();
            Ray rayOuter(PointOnScreen, rayDirectionOuter);

             value +=  triangle::getColourSoft(rayOuter, rayInner, root, vertices, faceVertices, faceTextures, textures,PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, innerLight.get_normal(), LightDirInner, &intersectionColours,  z, innerLight, outerLight);        
        }

        for(int z = 0; z<testIterations; z++){
            PixelColourSum += intersectionColours[z]; 
        }
        for(int z = 0; z<testIterations; z++){ //determine if one ray differs significantly
            if(((intersectionColours)[z]/PixelColourSum >0.05)&&(PixelColourSum>0)){
                adaptive = true;
            }
        } 

        if(adaptive==true){
            #pragma omp parallel for 
            for (int l=0; l<iterations;l++){
                vector3 PointOnInnerLight = innerLight.PointOnSource();
                vector3 rayDirection(PointOnInnerLight.x()-PointOnScreen.x(), PointOnInnerLight.y()-PointOnScreen.y(), PointOnInnerLight.z()-PointOnScreen.z());
                vector3 LightDirInner(PointOnScreen.x() - PointOnInnerLight.x(), PointOnScreen.y() - PointOnInnerLight.y(), PointOnScreen.z()-PointOnInnerLight.z());
                LightDirInner.normalize();
                rayDirection.normalize();
                Ray rayInner(PointOnScreen, rayDirection);

                vector3 PointOnOuterLight =outerLight.PointOnSource();
                vector3 rayDirectionOuter(PointOnOuterLight.x()-PointOnScreen.x(), PointOnOuterLight.y()-PointOnScreen.y(), PointOnOuterLight.z()-PointOnScreen.z());
                vector3 LightDirOuter(PointOnScreen.x() - PointOnOuterLight.x(), PointOnScreen.y() - PointOnOuterLight.y(), PointOnScreen.z()-PointOnOuterLight.z());
                LightDirOuter.normalize();
                rayDirectionOuter.normalize();
                Ray rayOuter(PointOnScreen, rayDirectionOuter);

                value +=  triangle::getColourSoft(rayOuter, rayInner, root, vertices, faceVertices, faceTextures, textures,PuppetTexture, PuppetTextureWidth, PuppetTextureHeight, innerLight.get_normal(), LightDirInner, &intersectionColours,  0, innerLight, outerLight);        
            }
        }
        delete [] intersectionColours; 

        //Using Monte Carlo, average values.
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
    std::ofstream image("puppet.bmp", std::ios::out| std::ios::binary); //write to bmp
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

	ObjFile::cleanUp(vertices, normals, textures, faceVertices, faceNormals, faceTextures);
    binarySearchTree::deleteTree(root);
    delete[]PuppetTexture;
    delete[] ScreenData;
    delete[] img;

    return 0;
}