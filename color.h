#ifndef COLOR_H
#define COLOR_H
using namespace std;

struct color {
	float r;
	float g;
	float b;
};

void col_copy (color *dest, color *src);

void col_init_fl (color *c, float r, float g, float b);

void col_init_uc (color *c, unsigned char r, unsigned char g, unsigned char b);

void col_add (color *res, color *a, color *b);

double col_diff (color a, color b);

void col_mul (color* res, color *a, double b);

void col_to_uchar_array (color *a, unsigned char *b);

void col_clamp (color *a);
#endif
