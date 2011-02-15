#include "stereoAC.h"
#include "libIO/io_tiff.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>

/// Ignore contrast, that is remove first PC.
static const bool bIgnoreContrast=true;

static bool loadImage(const char* name, LWImage<float>& im) {
    size_t nx, ny;
    im.data = read_tiff_f32_gray(name, &nx, &ny);
    im.w = nx; im.h = ny;
    if(! im.data)
        std::cerr << "Unable to load image file " << name << std::endl;
    return (im.data!=0);
}

/// Matrix mapping patch to principal components
static LWImage<float> principalComponents() {
#include "dataStereo/pca_basis.dat"
    int psize = sizeof(pca)/sizeof(*pca);
    int npcskept = sizeof(pca[0])/sizeof(*pca[0]);
    if(bIgnoreContrast)
        --npcskept;
    LWImage<float> PCs = alloc_image<float>(npcskept, psize);
    float* pc = PCs.data;
    for(int i=0; i < psize; i++)
        for(int j=0; j<npcskept; j++)
            *pc++ = pca[i][j+(bIgnoreContrast? 1: 0)];
    return PCs;
}

int main (int argc, char** argv)
{
    if(argc != 6 && argc != 7) {
        std::cerr << "Usage: " << argv[0] << " imgIn1 imgIn2 dMin dMax"
                  << " dispMapInc [dispMapMax]" << std::endl;
        return 1;
    }

    // Read input
    LWImage<float> im1(0,0,0), im2(0,0,0);
    if(! (loadImage(argv[1],im1) && loadImage(argv[2],im2)))
        return 1;

     // Read dispmin dispmax
    int dMin, dMax;
    std::istringstream f(argv[3]), g(argv[4]);
    if(! ((f>>dMin).eof() && (g>>dMax).eof())) {
        std::cerr << "Error reading dMin or dMax" << std::endl;
        return 1;
    }

    std::cout << "Meaningful Matches..." <<std::endl;
    /*-- VARIABLES --*/
    float epsNFA = 1.;
    epsNFA = static_cast<float>(im1.h);

    LWImage<float> PCs = principalComponents();
    float* dispInc = new float[im1.w*im1.h];
    float* dispMax = (argc>6)? new float[im1.w*im1.h]: 0;
    stereoAC(im1, im2, PCs, dMin, dMax, epsNFA, dispInc, dispMax);
    if(write_tiff_f32(argv[5], dispInc, im1.w, im1.h, 1) != 0) {
        std::cerr << "Error writing file " << argv[5] << std::endl;
        return 1;
    }
    if(dispMax != 0 && write_tiff_f32(argv[6], dispMax, im1.w, im1.h, 1) != 0) {
        std::cerr << "Error writing file " << argv[6] << std::endl;
        return 1;
    }

    free(PCs.data);
    free(im1.data);
	free(im2.data);
	delete [] dispInc;
    delete [] dispMax;

    return 0;
}
