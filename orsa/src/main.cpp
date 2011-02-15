#include "orsa.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#define _USE_MATH_DEFINES // For Windows
#include <cmath>
#include <ctime>

// Lexicographical ordering of matches. Used to remove duplicates.
static bool operator<(const Match& m1, const Match& m2)
{
    if(m1.x1 < m2.x1) return true;
    if(m1.x1 > m2.x1) return false;

    if(m1.y1 < m2.y1) return true;
    if(m1.y1 > m2.y1) return false;

    if(m1.x2 < m2.x2) return true;
    if(m1.x2 > m2.x2) return false;

    return (m1.y2 < m2.y2);
}

static bool operator==(const Match& m1, const Match& m2)
{
    return (m1.x1==m2.x1 && m1.y1==m2.y1 &&
            m1.x2==m2.x2 && m1.y2==m2.y2);
}

int main(int argc, char** argv)
{
    if(argc != 10) {
        std::cerr << "Usage: " << argv[0] << " w h match.txt good_match.txt ntrials verb noseed mode stop" <<std::endl;
        std::cerr << "w: width of image" <<std::endl;
        std::cerr << "h: height of image" <<std::endl;
        std::cerr << "match.txt: x1 y1 x2 y2 for each line" <<std::endl;
        std::cerr << "good_match.txt: good matchings (x1 y1 x2 y2 for each line)" <<std::endl;
        std::cerr << "ntrials: maximum number of ransac trials" <<std::endl;
        std::cerr << "verb: verbose mode (1 enabled, 0 disabled)" <<std::endl;
        std::cerr << "seed: random seed (0=reinitialize)" <<std::endl;
        std::cerr << "mode: 0=all 1=ransac 2=optimized ransac (ORSA) 3=automatic" <<std::endl;
        std::cerr << "stop: stop when first meaningful F is found (1 enabled, 0 disabled)" <<std::endl;
        return 1;
    }

    int width = 0, height = 0; // dimensions of image
    int ntrials = 0;           // maximum number of ransac trials
    bool verb = false;         // verbose
    unsigned long seed = 0;    // seed value (0=reinitialize)
    int mode = -1;    // 0=all 1=ransac 2=optimized ransac (ORSA) 3=automatic
    bool stop = false;         // stop when first meaningful F is found   

    if(! (std::istringstream(argv[1]) >> width).eof()) width = 0;
    if(! (std::istringstream(argv[2]) >> height).eof()) height = 0;
    if(width <=0 || height <= 0) {
        std::cerr << "Wrong dimensions of image" << std::endl;
        return 1;
    }

    std::vector<Match> match;
    if(! loadMatch(argv[3],match)) {
        std::cerr << "Failed reading " << argv[3] << std::endl;
        return 1;
    }

    if(! (std::istringstream(argv[5]) >> ntrials).eof() || ntrials <= 0) {
        std::cerr << "ntrials should be greater than 0" << std::endl;
        return 1;
    }

    if(! (std::istringstream(argv[6]) >> verb).eof()) {
        std::cerr << "verb can only be 0 or 1" << std::endl;
        return 1;
    }

    if(! (std::istringstream(argv[7]) >> seed).eof()) {
        std::cerr << "seed must be a non-negative integer value" << std::endl;
        return 1;
    }

    if(! (std::istringstream(argv[8]) >> mode).eof() || mode < 0 || mode > 3) {
        std::cerr << "mode can only be 0, 1, 2, or 3" << std::endl;
        return 1;
    }

    if(! (std::istringstream(argv[9]) >> stop).eof()) {
        std::cerr << "stop can only be 0 or 1" << std::endl;
        return 1;
    }

    // Initialize random seed if necessary
    if(seed == 0) {
        seed = (long int)time(NULL);
        if(verb)
            std::cout << "seed: " << seed << std::endl; // Useful for debugging
    }
    srand(seed);

    // Remove duplicates (frequent with SIFT)
    std::sort(match.begin(), match.end());
    std::vector<Match>::iterator end = std::unique(match.begin(), match.end());
    if(end != match.end()) {
        if(verb)
            std::cout << "Remove " << std::distance(end,match.end())
                      << "/" << match.size() << " duplicate matches"<<std::endl;
        match.erase(end, match.end());
    }

    // Normalize coordinates
    std::vector<Match> matchBackup(match);
    float nx = (float)width;
    float ny = (float)height;
    float norm = 1.0f/sqrt((float)(nx*ny));
    for(size_t i=0; i<match.size(); i++) {
        match[i].x1 =  (match[i].x1-0.5f*nx)*norm;
        match[i].y1 =  (match[i].y1-0.5f*ny)*norm;
        match[i].x2 =  (match[i].x2-0.5f*nx)*norm;
        match[i].y2 =  (match[i].y2-0.5f*ny)*norm;
    }
    libNumerics::matrix<float> N(3,3); // Normalization matrix
    N = 0;
    N(0,0) = N(1,1) = norm; N(2,2) = 1.0f;
    N(0,2) = -0.5f*nx*norm;
    N(1,2) = -0.5f*ny*norm;

    // log proba of a uniform point in image within a band of 1 pixel from line
    float logalpha0 = log10(2.0f)+0.5f*log10((nx*nx+ny*ny)/float(nx*ny));
    std::vector<size_t> inliers;
    float error;

	libNumerics::matrix<float> F = orsa(match, ntrials, verb, mode, stop, logalpha0, inliers, error);
    error /= norm;
    if(verb) {
        std::cout << "F= " << N.t()*F*N <<std::endl; // Denormalization
        std::cout << "Geometric error threshold: " << error <<std::endl;
    }

    // Write the good matchings into a file
    std::vector<Match> good_match;
    std::vector<size_t>::const_iterator it = inliers.begin();
    for(; it != inliers.end(); it++)
        good_match.push_back(matchBackup[*it]);

    if(! saveMatch(argv[4], good_match)) {
        std::cerr << "Failed saving good matchings into " <<argv[4] <<std::endl;
        return 1;
    }

    return 0;
}
