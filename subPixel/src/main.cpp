#include "subpixel.h"

#include "libIO/io_tiff.h"
#include "libLWImage/LWImage.h"

#include <iostream>
#include <cassert>

/*-- VARIABLES --*/
static const int NWIN = 4;

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
    if(argc != 5) {
        std::cerr << "Usage: " << argv[0] << " imgIn imgIn2 dispMapIn dispMapOut" << std::endl;
        return 1;
    }

    LWImage<float> im1(0,0,0), im2(0,0,0), disp(0,0,0);
    if(!(loadImage(argv[1],im1) &&
         loadImage(argv[2],im2) &&
         loadImage(argv[3],disp)))
        return 1;

    std::cout << "Subpixel refinement..." <<std::endl;

#include "dataStereo/prolate.dat"
    assert((4*NWIN+1)*(4*NWIN+1) == sizeof(prolate)/sizeof(*prolate));

    float* subDisp = new float[im1.w*im1.h];   

    refine_subpixel_accuracy(im1.data, im2.data, disp.data, subDisp,
                             prolate,4*NWIN+1, im1.w,im1.h, im2.w,im2.h);
    write_tiff_f32(argv[4], subDisp, im1.w, im1.h, 1);

    free(im1.data);
	free(im2.data);
	free(disp.data);
	delete [] subDisp;

    return 0;
}
