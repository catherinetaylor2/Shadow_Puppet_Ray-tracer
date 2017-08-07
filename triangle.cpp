#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "triangle.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include "search_tree.hpp"
#include "scene.hpp"

triangle::triangle(vector3 Vertex1, vector3 Vertex2, vector3 Vertex3){ //triangle constructor
    vertex1.setValue(Vertex1.x(), Vertex1.y(), Vertex1.z()); //triangle vertices
    vertex2.setValue(Vertex2.x(), Vertex2.y(), Vertex2.z());
    vertex3.setValue(Vertex3.x(), Vertex3.y(), Vertex3.z());

    vector3 Normal = vector3::crossproduct(vector3::subtract(vertex2, vertex1), vector3::subtract(vertex3, vertex1));
    Normal.normalize();
    normal.setValue(Normal.x(), Normal.y(), Normal.z()); 
    planeConstant = vector3::dotproduct(normal, vertex1); //n.x = d
}

float triangle::RayTriangleIntersection(Ray R){
    float r = vector3::dotproduct(normal, R.get_direction()); //check if intersects with plane
    if (fabs(r) < 0.000000001f){
              return 0;
    }
    float t = (planeConstant - vector3::dotproduct(normal,R.get_origin()))/r; //t where ray-plane intersection occurred
    vector3 intersectionPoint = vector3::add(R.get_origin(), vector3::ScalarMultiply(t,  R.get_direction())); //OuterPOI

    if( //test if inside triangle
    (vector3::dotproduct(vector3::crossproduct(vector3::add(vertex2, vector3::ScalarMultiply(-1, vertex1)), vector3::add(intersectionPoint, vector3::ScalarMultiply(-1, vertex1))), normal)>=-0.00001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::add(vertex3, vector3::ScalarMultiply(-1, vertex2)), vector3::add(intersectionPoint,vector3::ScalarMultiply(-1, vertex2))), normal)>=-0.00001f)&&
    (vector3::dotproduct(vector3::crossproduct(vector3::add(vertex1, vector3::ScalarMultiply(-1, vertex3)), vector3::add(intersectionPoint, vector3::ScalarMultiply(-1, vertex3))), normal)>=-0.00001f))
    {
        return t;
    }
    return 0;
}
 float triangle::getPOI(search_tree* root, float*vertices, Ray R, int* faces, int*minValue, int**k){ //use bounding boxes to find intersection.
    Bounding_box B_root(root->parameters[0],root->parameters[1], root->parameters[2],root->parameters[3],root->parameters[4],root->parameters[5]);
    std::vector<int> output; //stores possible triangles.
    float t_min = infinity;
	int index1, index2, index3;
    output.clear();
    if(B_root.ray_box_intersection(R)==1){
        search_tree::traverse_tree(root, R, &output); //test each bounding box.
    }
    *k = new int[output.size()+1];
    (*k)[0] = -1;
    if (output.size()>1){
        (*k)[0]=output.size();
        for(int g=1; g<(*k)[0]+1;g++){
            (*k)[g] = output[g-1]; //fill output with possible trianglws
        }
    }
    if( ((*k)[0]!=-1)&&((*k)[0]>0)){ //test each triangle in output/
        float* t_values = new float[(*k)[0]], t;
        int index;
        for (int z=1; z<(*k)[0]+1; z++){
            index = (*k)[z];
            index1 = faces[3*index] -1, index2 = faces[3*index+1]-1, index3 = faces[3*index+2] -1 ;
            vector3 Vertex1(vertices[3*index1], vertices[3*index1+1], vertices[3*index1+2]);
            vector3 Vertex2(vertices[3*index2], vertices[3*index2+1], vertices[3*index2+2]);
            vector3 Vertex3(vertices[3*index3], vertices[3*index3+1], vertices[3*index3+2]);
            triangle tri(Vertex1, Vertex2, Vertex3);
            t = tri.RayTriangleIntersection(R);
            t_values[z-1]=t;
        }
        for (int z=0; z<(*k)[0]; z++){ //find closest intersection.
            if ((t_values[z]>0)){
                vector3 xyz = vector3::add( R.get_origin(), vector3::ScalarMultiply(t_values[z],R.get_direction()));
                if (t_values[z]<t_min){
                    t_min =t_values[z];
                    *minValue = z;
                }
            }
        }
        delete t_values;
    }
    return t_min; //value of t when intersection occurs
 }

void triangle::getTextureValue(int triangleIndex, int* faceVertices, float *V, Ray R, unsigned char* puppetTexture, int* faceTextures, float* Textures, int puppetWidth, int puppetHeight, float **colour)
{ //use cramers rule to find baryentric coords.
    float denominator;    
    int index1, index2, index3;

    index1 = faceVertices[3*triangleIndex] -1, index2 = faceVertices[3*triangleIndex+1]-1, index3 = faceVertices[3*triangleIndex+2] -1 ; 
    vector3 point1(V[3*index1], V[3*index1+1], V[3*index1+2]);
    vector3 point2(V[3*index2], V[3*index2+1], V[3*index2+2]);
    vector3 point3(V[3*index3], V[3*index3+1], V[3*index3+2]);
    vector3 T = vector3::add(R.get_origin(), vector3::ScalarMultiply(-1, point1));
    vector3 E1 = vector3::add(point2, vector3::ScalarMultiply(-1, point1));
    vector3 E2 = vector3::add(point3, vector3::ScalarMultiply(-1, point1));
    denominator = vector3::dotproduct( vector3::crossproduct(R.get_direction(), E2), E1);
    vector3 M (vector3::dotproduct(vector3::crossproduct(T, E1), E2),vector3::dotproduct(vector3::crossproduct(R.get_direction(), E2), T),vector3::dotproduct( vector3::crossproduct(T, E1), R.get_direction()));
    vector3 tuv = vector3::ScalarMultiply(1.0f/(float)denominator, M);
    
    float *barycentric = new float[3]; //store barycentric coords
    (barycentric)[0] = 1.0f-(tuv.y()+tuv.z());
    (barycentric)[1] = tuv.y();
    (barycentric)[2] = tuv.z();

    int t_index1 = faceTextures[3*triangleIndex]-1, t_index2 = faceTextures[3*triangleIndex+1]-1, t_index3 = faceTextures[3*triangleIndex+2]-1; //get texture values from obj
    float vt_1x = Textures[2*t_index1], vt_1y = Textures[2*t_index1+1], vt_2x = Textures[2*t_index2], vt_2y = Textures[2*t_index2+1], vt_3x = Textures[2*t_index3], vt_3y = Textures[2*t_index3+1];

    float u_coord, v_coord, alpha, beta, Vertex12r, Vertex12g, Vertex12b, Vertex34r, Vertex34g, Vertex34b;
    int Vertex1x,Vertex1y, Vertex2x, v4y;
    u_coord = (barycentric[0]*vt_1x +barycentric[1]*vt_2x+barycentric[2]*vt_3x)*puppetWidth; //use barycentric coords to find position inside triangle
    v_coord = (barycentric[0]*vt_1y +barycentric[1]*vt_2y+barycentric[2]*vt_3y)*puppetHeight;
    Vertex1x = (int)floor(u_coord); //find corner pixel values
    Vertex1y = (int)ceil(v_coord);
    Vertex2x = (int)ceil(u_coord);
    v4y = (int)floor(v_coord);
    if (Vertex1x<0){
        Vertex1x=0;
    }
    if (Vertex2x<0){
        Vertex2x=0;
    }
    if (Vertex1y<0){
        Vertex1y=0;
    }
    if (v4y<0){
        v4y=0;
    }

    alpha = (float)(u_coord - (Vertex2x - Vertex1x)*Vertex1x)/(float) (Vertex2x - Vertex1x);
    beta = (float)(v_coord - (Vertex1y - v4y)*v4y)/(float) (Vertex1y - v4y);

    if (alpha >1){
        alpha=1;
    }
    if(beta>1){
        beta =1;
    }
//bilinear interpolation to find texture value
    Vertex12r = (1-alpha)*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex1x] + alpha*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex2x]; 
    Vertex12g = (1-alpha)*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex1x+1] + alpha*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex2x+1];
    Vertex12b = (1-alpha)*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex1x+2] + alpha*puppetTexture[Vertex1y*puppetWidth*3 + 3*Vertex2x+2];    
    Vertex34r = (1-alpha)*puppetTexture[v4y*puppetWidth*3 + 3*Vertex1x] + alpha*puppetTexture[v4y*puppetWidth*3 + 3*Vertex2x];
    Vertex34g = (1-alpha)*puppetTexture[v4y*puppetWidth*3 + 3*Vertex1x+1] + alpha*puppetTexture[v4y*puppetWidth*3 + 3*Vertex2x+1];
    Vertex34b = (1-alpha)*puppetTexture[v4y*puppetWidth*3 + 3*Vertex1x+2] + alpha*puppetTexture[v4y*puppetWidth*3 + 3*Vertex2x+2];

    (*colour)[0] = (1-beta)*Vertex12r + beta*Vertex34r; //colour of texture point
    (*colour)[1] = (1-beta)*Vertex12g + beta*Vertex34g;
    (*colour)[2] = (1-beta)*Vertex12b + beta*Vertex34b;

  delete[] barycentric;
}
float triangle::getColour(Ray R, search_tree* root, float*vertices, int* faceVertices, int*faceTextures, float*Textures,  unsigned char* puppetTexture, int puppetWidth, int puppetHeight, vector3 planeNormal, vector3 LightDirection, float**colours, int index){
    int minValue = -1, *k;
    float value;
    float t_min = triangle::getPOI(root, vertices, R,faceVertices, &minValue, &k); //test for intersection with quad
    if(minValue!=-1){ //if intersects with quad
        int triangle = k[minValue+1]; //triangle which has intersected with ray.
        float* colour = new float[3];
        triangle::getTextureValue(triangle, faceVertices, vertices, R, puppetTexture, faceTextures, Textures, puppetWidth, puppetHeight, &colour); //find value of texture at OuterPOI
        vector3 OuterPOI = vector3::add(R.get_origin(), vector3::ScalarMultiply(t_min,  R.get_direction()));  
        vector3 OuterDistVec = vector3::add(OuterPOI, vector3::ScalarMultiply(-1, R.get_origin()));
        float OuterDist = sqrt(vector3::dotproduct(OuterDistVec,OuterDistVec));
        float alpha = fabs(OuterDist)/75.0f; //OuterDistance function for level of blending

        if((colour[0]<10)&&(colour[1]<10)&&(colour[2]<10)){ //if intersects with puppet
            #pragma omp critical
            value = std::min(alpha,1.3f*pow(vector3::dotproduct(planeNormal, LightDirection),50.0f)+0.4f); //clamps shadow value at screen colour
        }
        else{ //intesects with quad but not puppet
            #pragma omp critical
            value = 1.3*pow(vector3::dotproduct(planeNormal, LightDirection),50.0f)+0.4f; //lighting model for background
        }  
        delete[] colour;
    }

    else{ //no intersections
        #pragma omp critical
        value = 1.3*pow(vector3::dotproduct(planeNormal, LightDirection),50.0f)+0.4f;              
    }
    (*colours)[index]=value; //save values for adaptive sampling
    delete[] k;
    return value;
    }
        
float triangle::getColourSoft(Ray rayOuter, Ray rayInner, search_tree* root, float* vertices, int* faceVertices, int*faceTextures, float*Textures, unsigned char* puppetTexture, int puppetWidth, int puppetHeight, vector3 planeNormal, vector3 L, float**colours, int index, light innerLight, light OuterLight){
    int minValueOuter = -1, *kInner, minValueInner = -1, *kOuter;
    float value;
    float t_minOuter = triangle::getPOI(root, vertices, rayOuter,faceVertices, &minValueOuter, &kInner);
    vector3 OuterPOI = vector3::add(rayOuter.get_origin(), vector3::ScalarMultiply(t_minOuter,  rayOuter.get_direction()));  
    vector3 OuterDistVec = vector3::add(OuterPOI, vector3::ScalarMultiply(-1, rayOuter.get_origin()));
    float OuterDist = sqrt(vector3::dotproduct(OuterDistVec,OuterDistVec));
    float alpha = fabs(OuterDist)/15.0f; //OuterDistance function for level of blending
                   
    float t_minInner = triangle::getPOI(root, vertices, rayInner,faceVertices, &minValueInner, &kOuter);     
    vector3 InnerOuterPOI = vector3::add(rayInner.get_origin(), vector3::ScalarMultiply(t_minInner,  rayInner.get_direction()));  
    vector3 InnerDistVec = vector3::add(InnerOuterPOI, vector3::ScalarMultiply(-1, rayInner.get_origin()));
    float InnerDist = sqrt(vector3::dotproduct(InnerDistVec,InnerDistVec));
    float InnerAlpha = fabs(InnerDist)/15.0f; //OuterDistance function for level of blending
   
    if(minValueInner!=-1){              
        if(minValueOuter!=-1){
            int triangle = kInner[minValueOuter+1];
            float* OuterColour = new float[3];
            triangle::getTextureValue(triangle, faceVertices, vertices, rayOuter, puppetTexture, faceTextures, Textures, puppetWidth, puppetHeight, &OuterColour); 
            int triangle2 = kOuter[minValueInner+1];
            float* InnerColour = new float[3];
            triangle::getTextureValue(triangle2, faceVertices, vertices, rayInner, puppetTexture, faceTextures, Textures, puppetWidth, puppetHeight, &InnerColour); 

            if((OuterColour[0]<10)&&(OuterColour[1]<10)&&(OuterColour[2]<10)){
                #pragma omp critical
                value =std::min(0.3f*OuterLight.get_illumination()*alpha,1.3f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f));
            }
            else if((InnerColour[0]<10)&&(InnerColour[1]<10)&&(InnerColour[2]<10)){
                #pragma omp critical
                value = std::min(0.85f*innerLight.get_illumination()*InnerAlpha,1.3f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f));                        
            }
            else{
                #pragma omp critical
                value =1.5f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f);
            }  
            delete[] OuterColour;
            delete[] InnerColour;
            delete[] kInner;
        }
        else{
            int triangle = kOuter[minValueInner+1];
            float* OuterColour = new float[3];
            triangle::getTextureValue(triangle, faceVertices, vertices, rayInner, puppetTexture, faceTextures, Textures, puppetWidth, puppetHeight, &OuterColour); 
            if((OuterColour[0]<10)&&(OuterColour[1]<10)&&(OuterColour[2]<10)){
                #pragma omp critical
                value = std::min(0.85f*innerLight.get_illumination(),1.3f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f));                        
            }
            else{
                #pragma omp critical
                value = 1.5f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f);
            }  
            delete[] OuterColour;
            delete[] kOuter;
        }
    }
    else{
        #pragma omp critical
        value = 1.3f*pow(vector3::dotproduct(innerLight.get_normal(), L),10.0f);              
    }
    (*colours)[index]=value; 
    return value;
}
        