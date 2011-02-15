#include "selfSimilar.h"
#include "libStereo/patch.h"
#include "libIO/nan.h"
#include <cfloat>
#include <cassert>

/// Check whether inter-image match is significantly higher than intra-image.
/// Inspired by SIFT matching criterion.
static bool check_sift(const LWImage<float>& im1, const LWImage<float>& im2,
                       int x1,int x2, int y, int lookr, int win, float ratioMax)
{
    assert(win<=y && y+win<im1.h && y+win<im2.h);
    assert(win<=x1 && x1+win<im1.w);
    assert(win<=x2 && x2+win<im2.w);
	float distmin = FLT_MAX; // Intra-image best SSD
	for(int i=-lookr; i<=lookr; i++) {
		int pos = x1+i;
		// Do not compare with itself and patches in position +/-1 pixel
		if(i!=0 && i!=1 && i!=-1 && win<=pos && pos+win<im1.w) {
            float dif = cssd(im1,x1,y, im1,pos,y, win);
            if(dif < distmin) distmin = dif;
		}
	}

	float dist =  cssd(im1,x1,y, im2,x2,y, win); // Inter-image SSD
	return (dist < ratioMax*distmin);
}

/// Filter out ambiguous patch matches from disparity map.
/// \a ratioMax of best and second best patch matches to be considered unique.
void selfSimilarFilter(const LWImage<float>& im1, const LWImage<float>& im2,
                       int lookr, int win,
                       float* disparity, float ratioMax)
{
    for(int y=0; y<im1.h; y++)
        for(int x=0; x<im1.w; x++, disparity++)
            if(is_number(*disparity) &&
               ! check_sift(im1, im2, x, x+*disparity, y, lookr, win, ratioMax))
                *disparity = NaN;
}
