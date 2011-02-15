#include "libTransform/TransformSize.h"
#include <float.h>

/// Constructor.
///
/// Rectangle is empty.
TransformSize::TransformSize()
: x1(FLT_MAX), y1(FLT_MAX), x2(-FLT_MAX), y2(-FLT_MAX)
{}

/// Add image of rectangle (0,0,\a w,\a h) by \a map to the size.
void TransformSize::add(const libNumerics::Homography& map, int w, int h)
{
    double x, y;

    x=0.0;   y=0.0;   map(x,y);
    if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
    if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

    x=w-1.0; y=0.0;   map(x,y);
    if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
    if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

    x=w-1.0; y=h-1.0; map(x,y);
    if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
    if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;

    x=0.0;   y=h-1.0; map(x,y);
    if(x < x1) x1 = (float)x; if(x > x2) x2 = (float)x;
    if(y < y1) y1 = (float)y; if(y > y2) y2 = (float)y;
}
