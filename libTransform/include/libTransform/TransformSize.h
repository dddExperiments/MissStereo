#ifndef TRANSFORMSIZE_H
#define TRANSFORMSIZE_H

#include "libNumerics/homography.h"

/// Compute bounding box of transformed rectangle(s).
///
/// This is used to compute the size of an image transformed by an
/// \c Homography, or of a mosaic.
class TransformSize
{
public:
    TransformSize();

    void add(const libNumerics::Homography& map, int w, int h);
    int x() const;
    int y() const;
    int w() const;
    int h() const;
private:
    float x1, y1, x2, y2; ///< Rectangle coordinates
};

/// Left coordinate of bounding box.
inline int TransformSize::x() const
{ return (int)x1; }

/// Top coordinate of bounding box.
inline int TransformSize::y() const
{ return (int)y1; }

/// Width of bounding box.
inline int TransformSize::w() const
{ return (x2 < x1)? 0: int(x2-x1+1); }

/// Height of bounding box.
inline int TransformSize::h() const
{ return (y2 < y1)? 0: int(y2-y1+1); }

#endif
