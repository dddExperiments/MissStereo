#include "libIO/io_png.h"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " image" << std::endl;
        return 1;
    }

    size_t w,h,c;
    unsigned char* data = read_png_u8(argv[1], &w,&h,&c);

    if(! data) {
        std::cerr << "Error reading image " << argv[1] << std::endl;
        return 1;
    }
    std::cout << w << " " << h << std::endl;

    return 0;
}
