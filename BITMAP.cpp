//Adapted from SOURCE: https://stackoverflow.com/questions/5751749/how-can-i-read-bmp-pixel-values-into-an-array
//

#include "BITMAP.hpp"
#include <iostream>
#include <fstream>

void fillBMPHeaders(BITMAP_File_Header *fileHeader, BITMAP_Info_Header *infoHeader, int width, int height){
    fileHeader->bfType = 0x4d42;
    fileHeader->bfSize =  54 + width*height*24;
    fileHeader->bfReserved1 = 0;
    fileHeader->bfReserved2 = 0;
    fileHeader->bfOffBits =54;

    infoHeader->biSize = 40;
    infoHeader->biWidth = width;
    infoHeader->biHeight =height ;
    infoHeader->biPlanes = 1;
    infoHeader->biBitCount = 24;
    infoHeader->biCompression = 0;
    infoHeader->biSizeImage = width*height*24;
    infoHeader->biXPelsPerMeter= 2835;
    infoHeader->biYPelsPerMeter= 2835;
    infoHeader->biClrUsed = 0;
    infoHeader->biClrImportant = 0;
}

void writeBMP(BITMAP_File_Header* fileHeader, BITMAP_Info_Header* infoHeader, std::ofstream *image){
    image->write((const char*)(&fileHeader->bfType), sizeof(fileHeader->bfType) );
    image->write((const char*)(&fileHeader->bfSize), sizeof(fileHeader->bfSize) );
    image->write((const char*)(&fileHeader->bfReserved1), sizeof(fileHeader->bfReserved1) );
    image->write((const char*)(&fileHeader->bfReserved2), sizeof(fileHeader->bfReserved2) );
    image->write((const char*)(&fileHeader->bfOffBits), sizeof(fileHeader->bfOffBits) );
    image->write((const char *)(infoHeader), sizeof(BITMAP_Info_Header));
}

unsigned char* readBMP(char* filename, int* imageWidth, int* imageHeight){
    int i;
    FILE* f = fopen(filename, "rb");
    if(f==nullptr){ //test if file exists
        return 0;
    }
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
    // extract image height and width from header
    *imageWidth = *(int*)&info[18];
    *imageHeight =abs(*(int*)&info[22]);
    int size = 3 * (*imageWidth) * (*imageHeight);
    unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
    fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
    fclose(f);
    for(i = 0; i < size; i += 3){
        unsigned char tmp = data[i];
        data[i] = data[i+2];
        data[i+2] = tmp;
    }
    return data;
}