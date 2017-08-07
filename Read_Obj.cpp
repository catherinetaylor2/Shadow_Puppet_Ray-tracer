#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include"Read_Obj.hpp"
#include "vec3.hpp"

#define infinity FLT_MAX

ObjFile::ObjFile(std::string name){ //constructor
	fn = name;
  if ((fopen(fn.c_str(), "r"))==nullptr){ //check file exists
    exist = false;
  }
  else{
    exist = true;
  }
}

void ObjFile::getVertices(float** vertices){ //find the vertices and store in array V
  char str[1000];
  float f1, f2, f3;
  FILE * myObject;
  myObject = fopen(fn.c_str(), "r"); //opens object file

  while (std::string(str) != "v"){ //scans file for specific arrangement of floats.
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
  do{  
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    NumberOfVertices += 1;
  }while (std::string(str) != "vt");

  *vertices= new float[3*NumberOfVertices];
  rewind(myObject); //go back to start

  while (std::string(str) != "v"){
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
  for(int i=0; i<3*NumberOfVertices; i+=3){ //save values in array.
    (*vertices)[i] = f1;
    (*vertices)[i+1] = f2;
    (*vertices)[i+2] = f3; 
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
  fclose(myObject);
}

void ObjFile::getNormals(float** N){ //as above with normals.
  char str[1000];
  float f1, f2, f3;
  std::string s = "a";
  FILE * myObject;
  int t;

  myObject = fopen(fn.c_str(), "r");
  while (s != "vn"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  do{  
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
    NumberOfNormals +=1;
  }while (s == "vn");

  rewind(myObject);
  *N = new float[3*NumberOfNormals];
  while (s != "vn"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i<3*NumberOfNormals; i+=3){
    (*N)[i] = f1;
    (*N)[i+1] = f2;
    (*N)[i+2]=f3;   
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
   fclose(myObject);
}

void ObjFile::getTextures(float ** VT){ //get texture values.
  char str[1000];
  float f1, f2, f3;
  std::string s = "a";
  FILE * myObject;
  int k_vt = 0;
  int t;
  myObject = fopen(fn.c_str(), "r");
  t = fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  s = str; 

  while (s != "v"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i< NumberOfVertices-1; i++){
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
  }
  do{  
    t=fscanf(myObject, "%s %f %f ", str, &f1, &f2);
    s = str;
    k_vt=k_vt+1;
  }while (s == "vt");

  fclose(myObject);
  myObject = fopen(fn.c_str(), "r");

  int number_of_text = k_vt-1;
  *VT = new float[2*number_of_text];
  while (s != "v"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i< NumberOfVertices-1; i++){
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
  }
  for(int i=0; i<2*number_of_text; i+=2){
    t=fscanf(myObject, "%s %f %f ", str, &f1, &f2);
    s=str;
    (*VT)[i] = f1;
    (*VT)[i+1]=f2;
  }
   fclose(myObject);
}

void ObjFile::getFaceData(int** faceVertices, int** faceNormals, int** faceTextures){ //find face vertices, normals and texture.
  char str[1000], c1, c2, c3, c4, c5, c6;
  float f1, f2, f3;
  int i1, i2, i3, i4, i5, i6, i7, i8, i9;
  std::string s = "a";
  FILE * myObject;
  int k_vf=0, t;

  myObject = fopen(fn.c_str(), "r");
  t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  s = str; 
  while (s != "vn"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i< NumberOfNormals-1; i++){
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
  }
  do{    
    s = str;
    k_vf=k_vf+1;
    t= fscanf(myObject, "%s %i %c %i %c %i %i %c %i %c %i %i %c %i %c %i", str, &i1, &c1, &i2, &c2, &i3, &i4, &c3, &i5, &c4, &i6, &i7, &c5, &i8, &c6, &i9);
  }while(t!=EOF);

  fclose(myObject);
  myObject = fopen(fn.c_str(), "r");
  t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  s = str; 
  while (s != "vn"){
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i< NumberOfNormals-1; i++){
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
  }
  NumberOfFaces = k_vf-1;
  std::cout<<"Number of mesh faces: "<<NumberOfFaces<<"\n";
  *faceVertices = new int[3*NumberOfFaces];
  *faceNormals = new int[3*NumberOfFaces];
  *faceTextures = new int[3*NumberOfFaces];

  for(int i=0; i<3*NumberOfFaces; i+=3){//store in arrays of integers
    t=fscanf(myObject, "%s %i %c %i %c %i %i %c %i %c %i %i %c %i %c %i", str, &i1, &c1, &i2, &c2, &i3, &i4, &c3, &i5, &c4, &i6, &i7, &c5, &i8, &c6, &i9);
    s=str;
    (*faceNormals)[i] = i3 -1;
    (*faceNormals)[i+1]=i6 -1;
    (*faceNormals)[i+2]=i9-1;
    (*faceVertices)[i]=i1-1;
    (*faceVertices)[i+1]= i4-1;
    (*faceVertices)[i+2]=i7-1;
    (*faceTextures)[i]=i2 -1;
    (*faceTextures)[i+1]=i5-1 ;
    (*faceTextures)[i+2]=i8 -1;
  }
   fclose(myObject);
}
void ObjFile::getMeshData(ObjFile mesh,int** faceVertices, int** faceNormals, int** faceTextures, float** texture_coords, float** normals, float** vertices, int* NumberOfFaces){
		mesh.getVertices(vertices);
		mesh.getTextures(texture_coords);
		mesh.getNormals(normals);
		mesh.getFaceData(faceVertices, faceNormals, faceTextures);
    NumberOfVertices = mesh.getNumberOfVertices();
	  *NumberOfFaces = mesh.getNumberOfFaces();
}
void ObjFile::cleanUp(float*vertices, float* normals, float* texture_coords,int* faceVertices, int* faceNormals, int* faceTextures){ //deletes pointers which are out of scope
  delete[] vertices;
  delete[] normals;
  delete[] texture_coords;
  delete[] faceVertices;
  delete[] faceNormals;
  delete[] faceTextures;
}
