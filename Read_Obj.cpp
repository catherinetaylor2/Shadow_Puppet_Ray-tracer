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
  if ((fopen(fn.c_str(), "r"))==nullptr){
    exist = false;
  }
  else{
    exist = true;
  }
}

void ObjFile::get_vertices(float** V){ //find the vertices and store in array V
  char str[1000];
  float f1, f2, f3;
  std::string s = "a";
  FILE * myObject;
  int k_v = 0, t;
  myObject = fopen(fn.c_str(), "r"); //opens object file

  while (s != "v"){ //scans file for specific arrangement of floats.
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  do{  
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
    k_v = k_v+1;
  }while (s != "vt");

  fclose(myObject);
  myObject = fopen(fn.c_str(), "r");
  number_of_vertices = k_v;
  *V = new float[3*number_of_vertices];
  while (s != "v"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i<3*number_of_vertices; i+=3){ //save values in array.
    (*V)[i] = f1;
    (*V)[i+1] = f2;
    (*V)[i+2]=f3; 
    t = fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
  fclose(myObject);
}

void ObjFile::get_normals(float** N){ //as above with normals.
  char str[1000];
  float f1, f2, f3;
  std::string s = "a";
  FILE * myObject;
  int k_vn=0, t;

  myObject = fopen(fn.c_str(), "r");
  while (s != "vn"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  do{  
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
    k_vn=k_vn+1;
  }while (s == "vn");

  fclose(myObject);
  myObject = fopen(fn.c_str(), "r");
  number_of_normals = k_vn;
  *N = new float[3*number_of_normals];
  while (s != "vn"){
    t= fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;  
  }
  for(int i=0; i<3*number_of_normals; i+=3){
    (*N)[i] = f1;
    (*N)[i+1] = f2;
    (*N)[i+2]=f3;   
    fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
  }
   fclose(myObject);
}

void ObjFile::get_texture(float ** VT){ //get texture values.
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
  for(int i=0; i< number_of_vertices-1; i++){
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
  for(int i=0; i< number_of_vertices-1; i++){
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

void ObjFile::get_face_data(int** face_vertex, int** face_normals, int** face_textures){ //find face vertices, normals and texture.
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
  for(int i=0; i< number_of_normals-1; i++){
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
  for(int i=0; i< number_of_normals-1; i++){
    t=fscanf(myObject, "%s %f %f %f" , str, &f1, &f2, &f3);
    s = str;
  }
  number_of_faces = k_vf-1;
  std::cout<<"Number of mesh faces: "<<number_of_faces<<"\n";
  *face_vertex = new int[3*number_of_faces];
  *face_normals = new int[3*number_of_faces];
  *face_textures = new int[3*number_of_faces];

  for(int i=0; i<3*number_of_faces; i+=3){//store in arrays of integers
    t=fscanf(myObject, "%s %i %c %i %c %i %i %c %i %c %i %i %c %i %c %i", str, &i1, &c1, &i2, &c2, &i3, &i4, &c3, &i5, &c4, &i6, &i7, &c5, &i8, &c6, &i9);
    s=str;
    (*face_normals)[i] = i3 -1;
    (*face_normals)[i+1]=i6 -1;
    (*face_normals)[i+2]=i9-1;
    (*face_vertex)[i]=i1-1;
    (*face_vertex)[i+1]= i4-1;
    (*face_vertex)[i+2]=i7-1;
    (*face_textures)[i]=i2 -1;
    (*face_textures)[i+1]=i5-1 ;
    (*face_textures)[i+2]=i8 -1;
  }
   fclose(myObject);
}
void ObjFile::get_mesh_data(ObjFile mesh,int** face_vertex, int** face_normals, int** face_textures, float** texture_coords, float** normals, float** vertices, int* number_of_faces){
		mesh.get_vertices(vertices);
		mesh.get_texture(texture_coords);
		mesh.get_normals(normals);
		mesh.get_face_data(face_vertex, face_normals, face_textures);
    number_of_vertices = mesh.get_number_of_vertices();
	  *number_of_faces = mesh.get_number_of_faces();
}
void ObjFile::clean_up(float*vertices, float* normals, float* texture_coords,int* face_vertex, int* face_normals, int* face_textures){ //deletes pointers which are out of scope
  delete[] vertices;
  delete[] normals;
  delete[] texture_coords;
  delete[] face_vertex;
  delete[] face_normals;
  delete[] face_textures;
}
