#include <color.h>
#include <image.h>
using namespace std;

void generate_headers (bmp_header *bmph, dib_header *dibh, int xsize, int ysize){
	int data_size = 3 * xsize * ysize;
	if ((xsize*3) % 4 != 0){
		data_size += ysize * (4-((xsize*3)%4));		// to account for padding
	}
	bmph->magic_number = (unsigned short) (66*256+77);//convert_ui((unsigned short) (66*256+77));
	bmph->size = (unsigned int) (54u + data_size);
	bmph->vacant = 0u;
	bmph->offset = 54u;
	
	dibh->header_bytes = 40u;
	dibh->bmp_width = (unsigned int) xsize;
	dibh->bmp_height = (unsigned int)ysize;
	dibh->col_planes = 1u;
	dibh->bpp = 24u;
	dibh->compression = 0u;
	dibh->data_size = (unsigned int)data_size;
	dibh->horiz_res = 2835u;
	dibh->vert_res = 2835u;
	dibh->palette_colors = 0u;
	dibh->imp_colors = 0u;
}

void insert_ushort (unsigned char *h, int place, unsigned short val){
	h[place] = val%256;
	h[place+1] = val/256;
}

void insert_uint (unsigned char *h, int place, unsigned int val){
	h[place] = val%256;
	val /= 256;
	h[place+1] = val%256;
	val /= 256;
	h[place+2] = val%256;
	val /= 256;
	h[place+3] = val;
}

void write_headers (bmp_header b, dib_header d, FILE *img){
	unsigned char h[54];
	h[0] = 66;
	h[1] = 77;
	insert_uint (h, 2, b.size);
	insert_uint(h, 6, b.vacant);
	insert_uint(h, 10, b.offset);

	insert_uint(h, 14, d.header_bytes);
	insert_uint(h, 18, d.bmp_width);
	insert_uint(h, 22, d.bmp_height);
	insert_ushort(h, 26, d.col_planes);
	insert_ushort(h, 28, d.bpp);
	insert_uint(h, 30, d.compression);
	insert_uint(h, 34, d.data_size);
	insert_uint(h, 38, d.horiz_res);
	insert_uint(h, 42, d.vert_res);
	insert_uint(h, 46, d.palette_colors);
	insert_uint(h, 50, d.imp_colors);

	fwrite(h, 1, 54, img);
}

void write_img (color *img, int xsize, int ysize, FILE *out){
	bmp_header b;
	dib_header d;
	generate_headers(&b,&d,xsize,ysize);
	write_headers(b,d,out);
	
	unsigned char col[3];
	for (int y=0; y<ysize; y++){
		for (int x=0; x<xsize; x++){
			col_to_uchar_array (&img[x*ysize+y], &col[0]);
			fwrite(&col[2], 1, 1, out);
			fwrite(&col[1], 1, 1, out);
			fwrite(&col[0], 1, 1, out);
		}
	}
	col[0] = (unsigned char) 0;
	if((xsize*3)%4 != 0){
		for (int i=0; i<(4-(xsize*3%4)); i++){
			fwrite(&col[0], 1, 1, out);
		}
	}
	fclose(out);
}
