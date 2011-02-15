#ifndef SELFSIMILAR_H
#define SELFSIMILAR_H

#include "libLWImage/LWImage.h"

void selfSimilarFilter(const LWImage<float>& im1, const LWImage<float>& im2,
                       int lookr, int win, float* disparity, float ratioMax=1);

#endif
