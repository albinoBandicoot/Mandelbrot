#include <color.h>
#include <cmath>
using namespace std;

void col_copy (color *dest, color *src){
	dest->r = src->r;
	dest->g = src->g;
	dest->b = src->b;
}

void col_init_fl (color *c, float r, float g, float b){
	c->r = r;
	c->g = g;
	c->b = b;
}

void col_init_uc (color *c, unsigned char r, unsigned char g, unsigned char b){
	c->r = r/255.0;
	c->g = g/255.0;
	c->b = b/255.0;
}

void col_add (color *res, color *a, color *b){
	res->r = a->r + b->r;
	res->g = a->g + b->g;
	res->b = a->b + b->b;
}

double col_diff (color a, color b){
	return (abs(a.r-b.r) + abs(a.g - b.g) + abs(a.b - b.b))/3.0;
}

void col_mul (color* res, color *a, double b){
	res->r = a->r*b;
	res->g = a->g*b;
	res->b = a->b*b;
}

void col_to_uchar_array (color *a, unsigned char *b){
	b[0] = (unsigned char) 255*a->r;
	b[1] = (unsigned char) 255*a->g;
	b[2] = (unsigned char) 255*a->b;
}

double max (double a, double b){
	if (a > b) return a;
	return b;
}

void col_clamp (color *a){
	a->r = max(0, a->r);
	a->g = max(0, a->g);
	a->b = max(0, a->b);
}
