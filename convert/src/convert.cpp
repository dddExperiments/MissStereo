#include "libIO/io_tiff.h"
#include "libIO/io_png.h"

#include <sstream>
#include <iostream>
#include <cfloat>
#include <cstdlib>

/// Put to true to have transparency layer
static const bool bTransparent=false;
/// Background color to put to invalid pixels
static const unsigned char bgColor[]={0,255,255};

// Find min/max
static void find_minmax(const float* im, int n, float& a, float& b)
{
    a=+FLT_MAX;
    b=-FLT_MAX;
    for(; n-- > 0; im++) {
        if(*im < a) a = *im;
        if(b < *im) b = *im;
    }
}

// Is there NaN in the image?
static bool has_nan(const float* im, int n)
{
    for(; n-- > 0; im++)
        if(*im != *im)
            return true;
    return false;
}

int main(int argc, char** argv)
{
    if(argc != 3 && argc != 5) {
        std::cerr << "Usage: " << argv[0] << " im_float.tif im.png [min max]" <<std::endl;
        return 1;
    }

    size_t w=0, h=0;
    float* im = read_tiff_f32_gray(argv[1], &w, &h);
    if(! im) {
        std::cerr << "Impossible to read float image " << argv[1] <<std::endl;
        return 1;
    }

    float a, b;
    if(argc > 4) {
        if(! (std::stringstream(argv[3]) >> a).eof()) {
            std::cerr << "Unable to interpret " << argv[3]
                      << " as min value"<<std::endl;
            return 1;
        }
        if(! (std::stringstream(argv[4]) >> b).eof()) {
            std::cerr << "Unable to interpret " << argv[4]
                      << " as max value"<<std::endl;
            return 1;
        }
    } else {
        find_minmax(im, w*h, a, b);
        std::cout << "min/max: " << a << " " << b << std::endl;
    }

    // Linear map from [a,b] to [256,0]: 256*(b-x)/(b-a)=(b-x)*(256/(b-a))
    if(b <= a) {
        a = -128.0f/a;
        b = 0;
    } else
        a = 256.0f / (b-a);

    const float* in = im;
    const int n = (bTransparent? 2: 3) * w*h;
    unsigned char* out = new unsigned char[n];
    unsigned char *red=out, *green=out+w*h, *blue=out+2*w*h;
    for(int i=w*h-1; i >= 0; i--, in++, red++) {
        float v = (b-*in)*a;
        if(v != v) v = bgColor[0];
        else if(v < 0.0f) v = 0.0f;
        else if(v > 255.0f) v = 255.0f;
        *red = static_cast<unsigned char>(v);
        if(bTransparent)
            *green++ = (*in==*in)? 255: 0;
        else {
            *green++ = (*in==*in)? *red: bgColor[1];
            *blue++  = (*in==*in)? *red: bgColor[2];
        }
    }

    bool bNaN = has_nan(im, w*h);
    const size_t channels = (bNaN? (bTransparent? 2: 3): 1);
    if(write_png_u8(argv[2], out, w, h, channels) != 0) {
        std::cerr << "Impossible to write png image " << argv[2] <<std::endl;
        return 1;
    }

    free(im);
    delete [] out;
    return 0;
}
