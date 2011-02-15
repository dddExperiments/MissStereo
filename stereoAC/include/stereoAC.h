#ifndef STEREOAC_H
#define STEREOAC_H

#include "libLWImage/LWImage.h"

void stereoAC(LWImage<float> im1, LWImage<float> im2, LWImage<float> PCs,
              int minDisp, int maxDisp, float epsNFA,
              float* dispInc, float* dispMax=0);

#endif
