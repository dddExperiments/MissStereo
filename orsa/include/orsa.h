#ifndef ORSA_H
#define ORSA_H

#include "libMatch/match.h"
#include "libNumerics/matrix.h"

libNumerics::matrix<float>
orsa(const std::vector<Match>& match,
     int t, bool verb, int mode, bool stop, float logalpha0,
     std::vector<size_t>& inliers, float& errorMax);

#endif
