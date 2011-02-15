#ifndef MAP_IMAGE_H
#define MAP_IMAGE_H

#include "libNumerics/homography.h"
#include "libLWImage/LWImage.h"
#include "libIO/nan.h"
#include <utility>

std::pair<int,int> map_image(LWImage<float> in,
                             libNumerics::Homography map,
                             LWImage<float>& im,
                             int order=1, bool adjustSize=true,
                             bool antiAlias=true,
                             float vOut=NaN);

#endif
