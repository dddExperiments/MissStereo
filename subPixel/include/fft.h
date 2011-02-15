#ifndef FFT_H
#define FFT_H

/*Add faster algorithm when size of signal is power of 2*/
void fft1d(float *Xr,float *Xi,float *Yr,float *Yi, int inverse, int n);
void fft2d(float *in_re, float *in_im,float *out_re, float *out_im, int i_flag, int width, int height);

void fft1dzoom(float *in, float *out, int simetry, int zoom, int n);

void fftzoom(float *in,float *out,float z, int width, int height);

#endif
