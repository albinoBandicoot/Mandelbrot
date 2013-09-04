#include <cmath>
#include <stdio.h>
#include <cmplx.h>
using namespace std;

Complex add (Complex a, Complex b){
	Complex res = {a.a+b.a, a.b+b.b};
	return res;
}

Complex mul (Complex a, Complex b){
	Complex res = {a.a*b.a - a.b*b.b, a.b*b.a + a.a*b.b};
	return res;
}

Complex mul_d (Complex a, double b){
	Complex res = {a.a*b, a.b*b};
	return res;
}

double dist (Complex x){
	return sqrt(x.a*x.a + x.b*x.b);
}

double dist_sqr (Complex x){
	return x.a*x.a + x.b*x.b;
}

void print (Complex a){
	printf("(%f", a.a);
	if (a.b < 0.0){
		printf(" - %fi)", abs(a.b));
	} else {
		printf(" + %f)i", a.b);
	}
}

/* 3-operand functions */

void add_3 (Complex *res, Complex *a, Complex *b){
	res->a = a->a + b->a;
	res->b = a->b + b->b;
}

void mul_3 (Complex *res, Complex *a, Complex *b){
	res->a = a->a*b->a - a->b*b->b;
	res->b = a->b*b->a + a->a*b->b;
}
