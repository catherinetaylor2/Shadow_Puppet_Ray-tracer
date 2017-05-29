//Adapted from SOURCE: https://stackoverflow.com/questions/5751749/how-can-i-read-bmp-pixel-values-into-an-array
//

#include "BITMAP.hpp"
#include <iostream>
#include <fstream>

void fill_bitmap_headers(BITMAP_File_Header *file_header, BITMAP_Info_Header *info_header, int width, int height){
    file_header->bfType = 0x4d42;
    file_header->bfSize =  54 + width*height*24;
    file_header->bfReserved1 = 0;
    file_header->bfReserved2 = 0;
    file_header->bfOffBits =54;

    info_header->biSize = 40;
    info_header->biWidth = width;
    info_header->biHeight =height ;
    info_header->biPlanes = 1;
    info_header->biBitCount = 24;
    info_header->biCompression = 0;
    info_header->biSizeImage = width*height*24;
    info_header->biXPelsPerMeter= 2835;
    info_header->biYPelsPerMeter= 2835;
    info_header->biClrUsed = 0;
    info_header->biClrImportant = 0;
}

void write_bitmap (BITMAP_File_Header* file_header, BITMAP_Info_Header* info_header, std::ofstream *image){
    image->write((const char*)(&file_header->bfType), sizeof(file_header->bfType) );
    image->write((const char*)(&file_header->bfSize), sizeof(file_header->bfSize) );
    image->write((const char*)(&file_header->bfReserved1), sizeof(file_header->bfReserved1) );
    image->write((const char*)(&file_header->bfReserved2), sizeof(file_header->bfReserved2) );
    image->write((const char*)(&file_header->bfOffBits), sizeof(file_header->bfOffBits) );
    image->write((const char *)(info_header), sizeof(BITMAP_Info_Header));
}

unsigned char* readBMP(char* filename, int* image_width, int* image_height){
    int i;
    FILE* f = fopen(filename, "rb");
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
    // extract image height and width from header
    *image_width = *(int*)&info[18];
    *image_height =abs(*(int*)&info[22]);
    int size = 3 * (*image_width) * (*image_height);
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