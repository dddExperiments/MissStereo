#ifndef _SPLINES_H_
#define _SPLINES_H_

#include "library.h"

float v(float *in,int x,int y,float bg, int width, int height);
void keys(float *c,float t,float a);
void spline3(float *c,float t);
void init_splinen(float *a,int n);
void splinen(float *c,float t,float *a,int n);

void finvspline(float *in,int order,float *out, int width, int height);

#endif
