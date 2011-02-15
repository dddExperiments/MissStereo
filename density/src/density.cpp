#include "libIO/io_tiff.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{
    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " im_float.tif" <<std::endl;
        return 1;
    }

    size_t w=0, h=0;
    float* im = read_tiff_f32_gray(argv[1], &w, &h);
    if(! im) {
        std::cerr << "Impossible to read float image " << argv[1] <<std::endl;
        return 1;
    }

    const float* in=im;
    int n=0;
    for(size_t i=w*h; i>0; i--, in++)
        if(*in == *in)
            ++n;
    std::cout << "Density: " << n << " /" << w*h << " = " << 100*n/(w*h) << "%"<<std::endl; 
    
    free(im);
    return 0;
}
