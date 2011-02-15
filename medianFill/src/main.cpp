#include "median_disparity.h"

#include "libIO/io_tiff.h"
#include "libLWImage/LWImage.h"

#include <iostream>

static bool loadImage(const char* name, LWImage<float>& im) {
    size_t nx, ny;
    im.data = read_tiff_f32_gray(name, &nx, &ny);
    im.w = nx; im.h = ny;
    if(! im.data)
        std::cerr << "Unable to load image file " << name << std::endl;
    return (im.data!=0);
}

int main (int argc, char** argv)
{
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " imgIn imgOut" << std::endl;
        return 1;
    }

    LWImage<float> im(0,0,0);
    if(! loadImage(argv[1],im))
        return 1;

    /*-------------------------------------------------------------*/
    /*-                  MEDIAN FILTER                            -*/
    /*-------------------------------------------------------------*/
    std::cout << "Median Filter..." <<std::endl;

    /*-- INITIALISATION --*/
    float* fill = new float[im.w*im.h];   

    median_disp(im.data, fill, 1, im.w, im.h);
    write_tiff_f32(argv[2], fill, im.w, im.h, 1);

    free(im.data);
    delete [] fill;

    return 0;
}
