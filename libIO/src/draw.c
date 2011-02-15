#include "libIO/draw.h"
#include <assert.h>

/** Draw a dashed horizontal line in image.
    \a data must have 3 non-interlaced channels.
    \a length is the length of 'on' sequence.
    \a gap is the length of 'off' sequence. */
void draw_horizontal_dashed_line(unsigned char* data, int w, int h,
                                 int y, int length, int gap, struct color c)
{
    unsigned char *r, *g, *b;
    int i, j;

    assert(0<=y && y<h);
    r = data + y*w + 0*w*h;
    g = data + y*w + 1*w*h;
    b = data + y*w + 2*w*h;

    for(i=0; i < w; i+=gap)
        for(j=0; j < length && i < w; j++, i++) {
            r[i] = c.r;
            g[i] = c.g;
            b[i] = c.b;
        }
}

/** Draw a cross in image.
    \a data must have 3 non-interlaced channels. */
void draw_cross(unsigned char* data, int w, int h,
                int x, int y, int halfLength, struct color c)
{
    unsigned char *r, *g, *b;
    int i;

    assert(0<=x && x<w && 0<=y && y<h);
    r = data + y*w+x + 0*w*h;
    g = data + y*w+x + 1*w*h;
    b = data + y*w+x + 2*w*h;

    for(i=-halfLength; i <= halfLength; i++) /* horizontal */
        if(0 <= x+i && x+i < w) {
            r[i] = c.r;
            g[i] = c.g;
            b[i] = c.b;
        }

    for(i=-halfLength; i <= halfLength; i++) /* vertical */
        if(0 <= y+i && y+i < h) {
            r[i*w] = c.r;
            g[i*w] = c.g;
            b[i*w] = c.b;
        }
}
