#ifndef MATCH_H
#define MATCH_H

#include <vector>

struct Match {
    Match() {}
    Match(float ix1, float iy1, float ix2, float iy2)
    : x1(ix1), y1(iy1), x2(ix2), y2(iy2) {}
    float x1, y1, x2, y2;
};

bool loadMatch(const char* nameFile, std::vector<Match>& match);
bool saveMatch(const char* nameFile, const std::vector<Match>& match);

#endif
