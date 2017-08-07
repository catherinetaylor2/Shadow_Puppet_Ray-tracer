//Adapted from SOURCE: https://stackoverflow.com/questions/5751749/how-can-i-read-bmp-pixel-values-into-an-array
//
#ifndef BITMAP_hpp
#define BITMAP_hpp

#include<iostream>
#include <fstream>

typedef int LONG;
typedef unsigned short WORD;
typedef unsigned int DWORD;

typedef struct BITMAP_FILE_HEADER {
  WORD  bfType;
  DWORD bfSize;
  WORD  bfReserved1;
  WORD  bfReserved2;
  DWORD bfOffBits;
} BITMAP_File_Header;

typedef struct BITMAP_INFO_HEADER {
  DWORD biSize;
  LONG  biWidth;
  LONG  biHeight;
  WORD  biPlanes;
  WORD  biBitCount;
  DWORD biCompression;
  DWORD biSizeImage;
  LONG  biXPelsPerMeter;
  LONG  biYPelsPerMeter;
  DWORD biClrUsed;
  DWORD biClrImportant;
} BITMAP_Info_Header;

void fillBMPHeaders(BITMAP_File_Header *fileHeader, BITMAP_Info_Header *infoHeader, int width, int height);
void writeBMP(BITMAP_File_Header *fileHeader, BITMAP_Info_Header *infoHeader, std::ofstream *image);
unsigned char* readBMP(char* filename, int* imageWidth, int* imageHeight); 

#endif