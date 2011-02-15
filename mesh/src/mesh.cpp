#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cfloat>
#include <cassert>

#include "libIO/io_tiff.h"
#include "libIO/io_png.h"

using namespace std;
using std::ofstream;

/// To use binary file format (more compact)
static const bool BINARY=false;

/// Uncomment next line to output a triangulation
//#define TRIANGULATION

/**this program implements a naive way of building a mesh from a disparity map:
the input is a float tiff image (a disparity map) and the output is a mesh in ascii ply file format
PLY file format can be read (among other) by paraview, meshlab. All modelers can also read PLY (blender, ogre...)
More information on PLY file format can be found at http://www.cc.gatech.edu/projects/large_models/ply.html
Usage:
disparitymap2mesh -i map.tiff -o output.ply
where map.tiff is the disparity map and output.ply is the mesh built from the map
*/

/**point structure containinng x,y,z and index of each pixel whose disparity has been computed*/
struct Point{
  int x;
  int y;
  float z;
  int index;
};

/**triangle structure containing the indices of the triangle vertices*/
struct Triangle{
  int i;
  int j;
  int k;
};

static std::string endian() {
    int i=0;
    *((unsigned char*)&i)=1;
    return (i==1? "little_endian": "big_endian");
}

int main(int argc, char **argv) {
     if(argc != 4) {
         cerr<<"Usage: " <<argv[0] <<" disp_float.tif text.png out.ply"<<endl;
         return EXIT_FAILURE;
     }

     //reading the image
      /*
     * number of columns, lines and channels of the image
     * size_t is the safe type,
     * because the image size is used as an array index
     */
    size_t nx, ny;

    /* read the image */
    float* img = read_tiff_f32_gray(argv[1], &nx, &ny);
    if (NULL == img) {
        cerr<<"failed to read the image "<<argv[1]<<endl;
        return EXIT_FAILURE;
    }
    
    size_t nx2, ny2;
    unsigned char* red = read_png_u8_rgb(argv[2], &nx2, &ny2);
    if(NULL == red) {
        cerr<<"failed to read the image "<<argv[2]<<endl;
        return EXIT_FAILURE;
    }
    unsigned char* green=red+nx2*ny2;
    unsigned char* blue=green+nx2*ny2;

    if(nx != nx2 || ny != ny2) {
        cerr<<"Error: images "<<argv[1]<< " and "<<argv[2]<< " must have same size"<<endl;
        return EXIT_FAILURE;
    }

    vector<Point> points;
    int index=0;
    float *px=img;
    float max=-FLT_MAX;
    float min=FLT_MAX;
    for(size_t y=0; y<ny; y++)
      for(size_t x=0; x<nx; x++, px++)
          if(*px==*px) {
              Point pt;
              pt.x=x;
              pt.y=y;
              pt.z=*px;
              pt.index=index++;
              points.push_back(pt);
              //record the highest and lowest disparity value to add a gray value to the point
              if(pt.z<min)
                  min=pt.z;
              if(pt.z>max)
                  max=pt.z;
          }	else {
              Point pt;
              pt.x=pt.y=pt.z=-1;
              pt.index=-1;
              points.push_back(pt);
          }

    free(img);
 
#ifdef TRIANGULATION
    vector<Triangle> triangles;
    for(size_t y=0; y+1<ny; y++) {
        vector<Point>::const_iterator it=points.begin()+y*nx;
        vector<Point>::const_iterator it2=it+nx;
        for(size_t x=0; x+1<nx; ++x, ++it, ++it2) {
            if(it->index!=-1 && (it+1)->index!=-1 && it2->index!=-1) {
                Triangle tr;
                tr.i=it->index;
                tr.j=(it+1)->index;
                tr.k=it2->index;
                triangles.push_back(tr);
            }
            if((it2+1)->index!=-1 && (it+1)->index!=-1 && it2->index!=-1) {
                Triangle tr;
                tr.i=(it2+1)->index;
                tr.j=(it+1)->index;
                tr.k=it2->index;
                triangles.push_back(tr);
            }
        }
    }
#endif

    //Let's initialize the ply file by writing the ply header
    ofstream out;
    out.open(argv[3]);
    out<<"ply"<<endl;
    if(BINARY)
        out<<"format binary_" << endian().c_str() << " 1.0"<<endl;
    else
        out<<"format ascii 1.0"<<endl;
    out<<"comment created by DisparityMap2Mesh"<<endl;
    out<<"element vertex "<<index<<endl;
    out<<"property float x"<<endl;
    out<<"property float y"<<endl;
    out<<"property float z"<<endl;
    out<<"property uchar red"<<endl;
    out<<"property uchar green"<<endl;
    out<<"property uchar blue"<<endl;
#ifdef TRIANGULATION
    out<<"element face "<<triangles.size()<<endl;
    out<<"property list uchar int vertex_index"<<endl;
#endif
    out<<"end_header"<<endl;

    //Now let's write all points in the file point is given a color corresponding to its elevation converted into a character value
    vector<Point>::const_iterator it=points.begin();
    const float d0=nx;
    const float deltaZ = nx/4.0f; // /10.0f;
    const float a=deltaZ/(1.0/d0 - 1.0/(max-min+d0));
    const float z0=deltaZ-a/d0;
    for(; it!=points.end(); ++it)
        if(it->index!=-1) {
            float z = a/(it->z-min+d0)+z0;
            size_t i=it->y*nx+it->x;
            if(BINARY) {
                float x=static_cast<float>(it->x);
                float y=static_cast<float>(ny-it->y-1);
                assert(sizeof(float)==4 && sizeof(unsigned char)==1);
                out.write((const char*)&x,4)
                   .write((const char*)&y,4)
                   .write((const char*)&z,4);
                out.write((const char*)&red[i],1)
                   .write((const char*)&green[i],1)
                   .write((const char*)&blue[i],1);
            } else
                out << it->x<<"\t"<<(ny-(size_t)it->y-1)<<"\t"<<z
                    <<"\t"<<(int)red[i]<<"\t"<<(int)green[i]<<"\t"<<(int)blue[i]<<endl;
        }

#ifdef TRIANGULATION
    //...and the triangles
    vector<Triangle>::const_iterator itt=triangles.begin();
    for(; itt!=triangles.end(); ++itt)
        if(BINARY) {
            assert(sizeof(int)==4);
            const unsigned char c=3;
            out.write((const char*)&c,1)
               .write((const char*)&itt->i,4)
               .write((const char*)&itt->j,4)
               .write((const char*)&itt->k,4);
        } else
        out<<3<<"\t"<<itt->i<<"\t"<<itt->j<<"\t"<<itt->k<<endl;
#endif

    out.close();
    return EXIT_SUCCESS;
}
