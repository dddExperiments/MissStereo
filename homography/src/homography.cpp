#include "libTransform/map_image.h"
#include "libNumerics/homography.h"
#include "libIO/io_png.h"
#include "libIO/io_tiff.h"
#include "libIO/nan.h"

#include <fstream>
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{
    if(argc != 4 && argc != 5) {
        std::cerr << "Usage: " <<argv[0] << " image_in H image_out [tiff32_out]"
                  << std::endl;
        return 1;
    }

    size_t nx, ny;
    LWImage<float> in(0,0,0);
    in.data = read_png_f32_rgb(argv[1], &nx, &ny);
    in.w = static_cast<int>(nx); in.h = static_cast<int>(ny); in.comps=3;
    if(! in.data) {
        std::cerr << "Error reading image " << argv[1] << std::endl;
        return 1;
    }

    libNumerics::Homography H;
    std::ifstream f(argv[2]);
    if((f >> H.mat()).fail()) {
        std::cerr << "Error reading homography file " << argv[2] << std::endl;
        return 1;
    }

    LWImage<float> out = alloc_image<float>(in.w, in.h, in.comps);
    map_image(in, H, out, 5, false, true);

    // Write only red channel in float for the moment
    if(argc >= 5 &&
       write_tiff_f32(argv[4], out.data, out.w, out.h, 1/*out.comps*/) != 0) {
        std::cerr << "Error writing file " << argv[4] << std::endl;
        return 1;
    }

    // Put in white invalid pixels
    for(int i=out.comps*out.w*out.h-1; i>=0; i--)
        if(! is_number(out.data[i]))
            out.data[i] = 255.0f;
    if(write_png_f32(argv[3], out.data, out.w, out.h, out.comps) != 0) {
        std::cerr << "Error writing file " << argv[3] << std::endl;
        return 1;
    }

    free(in.data);
    free(out.data);

    return 0;
}
