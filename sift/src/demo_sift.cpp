#include "demo_lib_sift.h"
#include "library.h"

#include "libIO/io_png.h"

#include <iostream>
#include <fstream>

bool load(char* nameFile, float*& pix, int& w, int& h)
{
    size_t sw,sh;
    pix = read_png_f32_gray(nameFile, &sw,&sh);
    if(! pix) {
        std::cerr << "Unable to load image file " << nameFile << std::endl;
        return false;
    }

	w = static_cast<int>(sw);
    h = static_cast<int>(sh);
    return true;
}

int main(int argc, char **argv)
{	
    if(argc != 4 && argc != 5) {
        std::cerr << "Usage: " << argv[0] << " imgIn imgIn2 fileOut [imgOut]"
                  << std::endl;
        return 1;
    }

	//////////////////////////////////////////////// Input
	int w1, h1;
    float* ipixels1;
    if(! load(argv[1], ipixels1, w1, h1))
        return 1;

	//////////////////////////////////////////////// Input
    int w2, h2;
    float* ipixels2;
    if(! load(argv[2], ipixels2, w2, h2))
        return 1;

	///////////////////////////////////////// Applying Sift
	siftPar siftparameters;
	default_sift_parameters(siftparameters);
    siftparameters.DoubleImSize=0;

	keypointslist keyp1, keyp2;
	compute_sift_keypoints(ipixels1,keyp1,w1,h1,siftparameters);
    std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints"<<std::endl;
	compute_sift_keypoints(ipixels2,keyp2,w2,h2,siftparameters);
    std::cout<< "sift:: 2nd image: " << keyp2.size() << " keypoints"<<std::endl;

	matchingslist matchings;
	compute_sift_matches(keyp1,keyp2,matchings,siftparameters);	
    std::cout << "sift:: matches: " << matchings.size() <<std::endl;

	//////////////////////////////////////////////////////////////// Save file with matches
    saveMatch(argv[3], matchings);

	//////////////////////////////////////////////// Output image containing line matches
    if(argc > 4) {
        int wo =  std::max(w1,w2);
        int ho = h1+h2;

        float *opixels = new float[wo*ho];
        for(int j = 0; j < h1; j++)
            for(int i = 0; i < w1; i++)  opixels[j*wo+i] = ipixels1[j*w1+i];
        for(int j = 0; j < h2; j++)
            for(int i = 0; i < w2; i++)  opixels[(h1 + j)*wo + i] = ipixels2[j*w2 + i];	

        //////////////////////////////////////////////////////////////////// Draw matches
        matchingslist::iterator ptr = matchings.begin();
        for(; ptr != matchings.end(); ++ptr)
            draw_line(opixels,
                      (int) ptr->x1, (int) ptr->y1,
                      (int) ptr->x2, (int) ptr->y2 + h1,
                      255.0f, wo, ho);

        ///////////////////////////////////////////////////////////////// Save imgOut	
        write_png_f32(argv[4], opixels, (size_t)wo, (size_t)ho, 1);
        delete[] opixels;
	}

	/////////////////////////////////////////////////////////////// Delete memory
	free(ipixels1);
	free(ipixels2);
    return 0;
}
