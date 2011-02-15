#include "libMatch/match.h"
#include "libNumerics/homography.h"
#include "libIO/draw.h"
#include "libIO/io_png.h"
#include "libLWImage/LWImage.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>

static const int NUM_LINES=15; ///< Number of lines to display in image
static const int LENGTH_ON=7; ///< Length of 'on' pattern in dashed line
static const int LENGTH_OFF=3; ///< Length of 'off' pattern in dashed line
static const color CYAN={0,255,255};
static const int CROSS_HALF_LENGTH = 3; ///< Half-length of cross for SIFT
static const color RED={255,0,0};

int main(int argc, char** argv)
{
    if(argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " in.png out.png match.txt left|right H.txt" <<std::endl;
        return 1;
    }

    // Read image
    size_t nx, ny;
    LWImage<unsigned char> in(0,0,0);
    in.data = read_png_u8_rgb(argv[1], &nx, &ny);
    in.w = static_cast<int>(nx); in.h = static_cast<int>(ny); in.comps=3;
    if(! in.data) {
        std::cerr << "Error reading image " << argv[1] << std::endl;
        return 1;
    }

    // Read correspondences
    std::vector<Match> match;
    if(! loadMatch(argv[3],match)) {
        std::cerr << "Failed reading " << argv[3] << std::endl;
        return 1;
    }

    // Read left/right image
    bool bLeft = (strcmp(argv[4],"left")==0);
    if(!bLeft && strcmp(argv[4], "right")!=0) {
        std::cerr << "Error: arg 4 must be 'left' or 'right'" << std::endl;
        return 1;
    }

    // Read homography
    libNumerics::Homography H;
    std::ifstream f(argv[5]);
    if((f >> H.mat()).fail()) {
        std::cerr << "Error reading homography file " << argv[2] << std::endl;
        return 1;
    }

    // Drawing lines
    int delta = in.h/(NUM_LINES+1);
    if(delta==0) delta=1;
    for(int i=delta; i < in.h; i+=delta)
        draw_horizontal_dashed_line(in.data, in.w, in.h, i,
                                    LENGTH_ON, LENGTH_OFF, CYAN);

    // Drawing SIFT
    std::vector<Match>::const_iterator it=match.begin();
    for(; it != match.end(); ++it) {
        double x=it->x1, y=it->y1;
        if(! bLeft) {
            x = it->x2; y=it->y2;
        }
        H(x,y);
        int ix = static_cast<int>(std::floor(x+0.5));
        int iy = static_cast<int>(std::floor(y+0.5));
        if(0<=ix && ix<in.w && 0<=iy && iy<in.h)
            draw_cross(in.data, in.w, in.h, ix, iy, CROSS_HALF_LENGTH, RED);
    }

    // Write image
    if(write_png_u8(argv[2], in.data, in.w, in.h, in.comps) != 0) {
        std::cerr << "Error writing file " << argv[3] << std::endl;
        return 1;
    }
    return 0;
}
