#include "libStereo\patch.h"

/// Sum of pixel values in patch centered on (i,j).
float sum(const LWImage<float>& im, int i, int j, int win)
{
    float s=0.0f;
	for(int y=-win; y<=win; y++) {
		const float* p = im.pixel(i-win, j+y);
        for(int x=-win; x<=win; x++)
            s += *p++;
    }
    return s;
}

/// Sum of Square Differences between patches centered on (i1,j1) and (i2,j2).
float ssd(const LWImage<float>& im1, int i1,int j1, 
          const LWImage<float>& im2, int i2,int j2, int win)
{
	float dist=0.0f;
	for(int j=-win; j<=win; j++) {
		const float* p1 = im1.pixel(i1-win, j1+j);
		const float* p2 = im2.pixel(i2-win, j2+j);
		for(int i=-win; i<=win; i++){
			float dif = (*p1++ - *p2++);
			dist += (dif*dif);
		}
	}
	return dist;
}

/// Centered Sum of Square Differences of patches of size 2*win+1.
float cssd(const LWImage<float>& im1,int i1,int j1,
           const LWImage<float>& im2,int i2,int j2, int win)
{
    float m = sum(im1,i1,j1,win)-sum(im2,i2,j2,win);
    int w = 2*win+1;
    return (ssd(im1,i1,j1, im2,i2,j2, win) - m*m/(w*w));
}
