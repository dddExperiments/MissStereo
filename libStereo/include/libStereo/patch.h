#ifndef PATCH_H
#define PATCH_H

#include "libLWImage/LWImage.h"

float  sum(const LWImage<float>& im, int i, int j,  int win);
float  ssd(const LWImage<float>& im1,int i1,int j1,
           const LWImage<float>& im2,int i2,int j2, int win);
float cssd(const LWImage<float>& im1,int i1,int j1,
           const LWImage<float>& im2,int i2,int j2, int win);

#endif
