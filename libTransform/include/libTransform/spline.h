#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "libLWImage/LWImage.h"
#include "libIO/nan.h"

bool prepare_spline(LWImage<float>& im, int order);
bool interpolate_spline(LWImage<float>& im, int order,
                        float x, float y,
                        float* out,
                        float paramKeys=-.5f);

#endif
