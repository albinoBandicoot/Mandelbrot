#ifndef COMPLEX_H
#define COMPLEX_H
#include <cmath>
#include <stdio.h>
using namespace std;

struct Complex {
	double a;
	double b;
};

Complex add (Complex a, Complex b);

Complex mul (Complex a, Complex b);

Complex mul_d (Complex a, double b);

double distance (Complex x);

double dist_sqr (Complex x);

void print (Complex a);

/* 3-operand functions */

void add_3 (Complex *res, Complex *a, Complex *b);

void mul_3 (Complex *res, Complex *a, Complex *b);

#endif
