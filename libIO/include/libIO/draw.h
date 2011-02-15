#ifndef DRAW_H
#define DRAW_H

#ifdef __cplusplus
extern "C" {
#endif

struct color {
    unsigned char r, g, b;
};

void draw_horizontal_dashed_line(unsigned char* data, int w, int h,
                                 int y, int length, int gap, struct color c);
void draw_cross(unsigned char* data, int w, int h,
                int x, int y, int halfLength, struct color c);

#ifdef __cplusplus
}
#endif

#endif
