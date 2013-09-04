//#ifndef IMAGE_H
//#define IMAGE_H

#include <color.h>
#include <stdio.h>
using namespace std;

typedef struct bmpheader {
	unsigned short magic_number;
	unsigned int size;
	unsigned int vacant;
	unsigned int offset;
} bmp_header;

typedef struct dibheader {
	unsigned int header_bytes;
	unsigned int bmp_width;
	unsigned int bmp_height;
	unsigned short col_planes;
	unsigned short bpp;
	unsigned int compression;
	unsigned int data_size;
	unsigned int horiz_res;
	unsigned int vert_res;
	unsigned int palette_colors;
	unsigned int imp_colors;
} dib_header;

void generate_headers (bmp_header *bmph, dib_header *dibh, int xsize, int ysize);

void insert_ushort (unsigned char *h, int place, unsigned short val);

void insert_uint (unsigned char *h, int place, unsigned int val);

void write_headers (bmp_header b, dib_header d, FILE *img);

void write_img (color *img, int xsize, int ysize, FILE *out);

//#endif
