#include "median_disparity.h"
#include "libIO/nan.h"

#include <algorithm>
#include <vector>
#include <stdlib.h>

void median_disp(float *disparity, float *median_disparity, int win, int nx, int ny)
{
    const int nb_pixels = nx*ny;
    const int ini = win;
    const int dim_patch = (2*win+1)*(2*win+1);
    const int fin_x = nx-ini;
    const int fin_y = ny-ini;

    const size_t pmin = (dim_patch+1)/2; //Half window
    int iter = 2;

    std::vector<float> window;

    for(int k=0; k<nb_pixels; k++)
        median_disparity[k] = disparity[k];

    while(iter--) {
        for(int y=ini; y<fin_y; y++){
            int ycolumn = y*nx;
            for(int x=ini; x<fin_x; x++){
                int POS = ycolumn + x;
                if( is_number(disparity[POS]) ) // Interpolation only
                    continue;

                window.clear();
                for(int j=-win; j<=win; j++) {
                    int columns = (y+j)*nx; 
                    for(int i=-win; i<=win; i++) {
                        int pos = columns + (x+i);
                        if( is_number(disparity[pos]) )
                            window.push_back(disparity[pos]);
                    }
                }

                size_t p = window.size();
                if(p >= pmin){
                    std::vector<float>::iterator it=window.begin()+(p-1)/2;
                    std::nth_element(window.begin(), it, window.end());
                    if((window.size()&1) == 0) // Even number
                        *it = (*it+*std::min_element(it+1, window.end()))/2;
                    median_disparity[POS] = *it;
                }
            }
        }
        for(int k=0; k<nb_pixels; k++) // Prepare for next iteration
            disparity[k] = median_disparity[k];
    }
}
