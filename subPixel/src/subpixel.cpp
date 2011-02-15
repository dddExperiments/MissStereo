#include "subpixel.h"
extern "C" {
#include "fft.h"
}
#include "libTransform/map_image.h"
#include "libIO/nan.h"
#include <stdlib.h>
#include <stdio.h>

/// Numerical accuracy, to avoid dividing by 0
static const float EPS=1.0e-5f;

/// Replace NaN pixels with value \a v.
static void correct_nan(float* data, int w, int h, float v=0.0f) {
    for(int i=w*h-1; i>=0; i--, data++)
        if(! is_number(*data))
            *data = v;
}

/*Euclidean weighted distance of patches*/
float window_kernel_distance(float *u0,float *u1,int i0,int j0,int i1,int j1,int width,int width2,float *kernel,int radius)
{
 
	int s,r,l;
	float *ptrk = &kernel[0];
	float dif , dist =0;   
	float *ptr0, *ptr1;

	for (s=-radius; s<= radius; s++){
	
		l = (j0+s)*width + (i0-radius);
		ptr0 = &u0[l];
	
		l = (j1+s)*width2 + (i1-radius);
		ptr1 = &u1[l];
	
		for(r=-radius;r<=radius;r++,ptr0++,ptr1++,ptrk++){
	  
			dif = (float) (*ptr0 - *ptr1);
			dist += (float) *ptrk*(dif*dif);		
		}
	}

	return dist;

}


/*Symmetrization of an image*/
void md_fim_sim(float *igray,float *ogray, int inverse,int width, int height)
{

  int x,y,size, nx,ny;
 
  nx = width;
  ny = height;
  size = nx*ny;
 
  if (inverse) {
    /* croping */
    if (nx & 1 || ny & 1) printf("WARNING: Non-even image dimensions\n"); 

    for (x=0; x < nx/2; x++)
      for (y= 0; y < ny/2; y++) 
	{
		ogray[y*(nx/2) + x] = igray[y*nx + x];
	}
     
  } else {

    /* symmetrization */
    for (x=nx;x--;)
      for (y=ny;y--;) 
	ogray[y*nx*2+x]
	  = ogray[y*nx*2+nx*2-1-x]
	  = ogray[(ny*2-1-y)*nx*2+x]
	  = ogray[(ny*2-1-y)*nx*2+nx*2-1-x]
	  = igray[y*nx+x];
  }

}


/*Zoom by simmetry and Fourier*/
void sim_fzoom(float *igray, float *ogray, int zoom, int width, int height)
{

	float *sigray;
	float *zsigray;
	sigray=(float *)malloc(4*width*height*sizeof(float));

	md_fim_sim(igray,sigray,0,width,height);

	zsigray=(float *)malloc(4*zoom*zoom*width*height*sizeof(float));
	fftzoom(sigray,zsigray,(float) zoom,2*width,2*height);

	md_fim_sim(zsigray, ogray, 1, 2*zoom*width,2*zoom*height);

	free(sigray);
	free(zsigray);
}


void fft1Dinterp_filtre(float *in,float *out,int z, int length)
{

  int x,adr, adrz;	

  int zlength = z * length;
  int hx = length / 2;


  /*filter in Fourier domain*/
  float filtre[16];
  float *aux_re, *aux_im;
  float *in_im;
  float *zaux_re, *zaux_im;

  filtre[0] = 0.104213f;   filtre[1] = 0.101741f;    filtre[2] = 0.0945092f;
  filtre[3] = 0.0830657f;  filtre[4] = 0.0682754f;   filtre[5] = 0.0512564f;
  filtre[6] = 0.0332955f;  filtre[7] = 0.0157505f;   filtre[8] = 0.0157505f;
  filtre[9] = 0.0332955f;  filtre[10] = 0.0512564f;  filtre[11] = 0.0682754f;
  filtre[12] = 0.0830657f; filtre[13] = 0.0945092f;  filtre[14] = 0.101741f;
  filtre[15] = 0.104213f;


  /*Auxiliar to compute fft of input signal*/
  aux_re=(float *)malloc(length*sizeof(float));
  aux_im=(float *)malloc(length*sizeof(float));

  /*Auxiliar as imaginary part of input signal*/
  in_im=(float *)malloc(length*sizeof(float));

  for(x=0;x < length; x++)
    in_im[x]=0.;

  /*Compute fft of input signal*/
  fft1d(in,in_im,aux_re,aux_im, 0, length);

  zaux_re=(float *)malloc(zlength*sizeof(float)); 
  zaux_im=(float *)malloc(zlength*sizeof(float));

  for( x=0; x < zlength; x++)
    {  zaux_re[x]=0.;  zaux_im[x]=0.;  }

  for (x = -hx; x <= hx;x++){

	adr  = (x+length )%length;
    	adrz = (x+zlength)%zlength;

    	zaux_re[adrz] = aux_re[adr]*filtre[adr];
    	zaux_im[adrz] = aux_im[adr]*filtre[adr];
  }

  fft1d(zaux_re,zaux_im,out,NULL,1, zlength);

 free(zaux_re);
 free(zaux_im);
 free(aux_re);
 free(aux_im);
 free(in_im);

}

void refine_subpixel_accuracy(float *input1,float *input2,float *disparity,float *subpixeldisparity, float* prolate, int prolwidth, int width, int height, int width2, int height2)
{

	/*Initialisation*/
	int x,i,y,k,l,x0,x1,pos_min;
	float cmin,y_p,y_q,y_r,x_min,dif;
	int zoom = 2; /*Initial zoom*/ 
	int dzoom = 32; /*Zoom for 1D interpolation of the quadratic distance*/ 

	int zwidth = zoom*width;    /*size images after zoom*/
	int zheight = zoom*height;
	int zwidth2 = zoom*width2;
	int zheight2 = zoom*height2;

	int kmin = -4;
	int kmax = 3;
 	int epilong = kmax-kmin+1;
	int sym_epilong = 2*epilong;
	int bound0 = dzoom*2;            /*Boundary for correlation*/
	int bound1 = dzoom*(epilong -1); /*Boundary for correlation*/
	float  Z = (float)(zoom*dzoom);
	float  T = (float) kmin / (float) zoom;
	float *correlation,*zcorrelation; 		
	float *zinput1;
	float *zinput2;


	int rad = (prolwidth - 1) / 2; /*radious of prolate patch.*/
	int rad_patch = rad /2; /*radious of the patch*/

	int pw0 = rad_patch - kmin/2;  /*Image boundary condition*/ 
 	int pw1 = width2 - pw0;             /*Image boundary condition*/

	for(i=0; i < width*height; i++) {subpixeldisparity[i] = NaN;}

	correlation = (float *)malloc(sym_epilong*sizeof(float)); /*2x samples for symmetrization*/	
	zcorrelation = (float *)malloc(sym_epilong*dzoom*sizeof(float));/*in order to interpolate disparity*/

    libNumerics::Homography H;
    H.setZoom(zoom,zoom);

	/*Initial zoom of 2 to avoid aliasing with symmetry*/
    correct_nan(input1, width,height);
	zinput1=(float *)malloc(zwidth*zheight*sizeof(float));
	//sim_fzoom(input1,zinput1,zoom, width,height);
    LWImage<float> out = make_image(zinput1,zwidth,zheight);
    map_image(make_image(input1,width,height), H, out, 5, false, false, 0.0f);

    correct_nan(input2, width2,height2);
	zinput2=(float *)malloc(zwidth2*zheight2*sizeof(float));
	//sim_fzoom(input2,zinput2,zoom, width2,height2);
    out = make_image(zinput2,zwidth2,zheight2);
    map_image(make_image(input2,width,height), H, out, 5, false, false, 0.0f);
	
	/*Refinate disparities*/

  	for(x = rad_patch+1; x < width -rad_patch -1; x++){
    		for( y= rad_patch + 1; y < height  - rad_patch -1; y++)

                if (is_number(disparity[y*width+x])){
      	
				/*Compute disparity*/
				l = y*width+x;
				x0 = x + disparity[l];  /*x0 is the corresponding pixel in the second image*/

				if (x0  >=  pw0 && x0 <= pw1 ){  /*boundary condition on the second image*/
				
					for(k=0; k < sym_epilong; k++) correlation[k] = 0.0;

					for(k = kmin ; k <= kmax; k++){

						x1 = 2*x0 + k;
		
						dif =  window_kernel_distance(zinput1,zinput2,2*x,2*y,x1,2*y,
								zwidth,zwidth2,prolate,rad);
					
						correlation[k - kmin] = dif; 
					
					}

					/*Symmetry*/
				 	for(k=epilong;k < sym_epilong; k++)
	 				  	correlation[k] = correlation[sym_epilong-k-1];			
			 		
					/*Interpolation of correlation*/
				 	fft1Dinterp_filtre(correlation,zcorrelation,dzoom,sym_epilong);

	  				/*Searching the min inside (bound0, bound1) i.e. positions +/-1 pixels*/
	  				cmin = zcorrelation[bound0];
	  				pos_min = bound0;

		  			for( k = bound0; k <= bound1; k++){	
							    			
						if (zcorrelation[k] < cmin) { cmin = zcorrelation[k]; pos_min =  k;}
					}
					
					/*Discard boundary*/
					if (pos_min != bound0 && pos_min != bound1) { 
					
						/*Final quadratic interpolation*/				
						y_p = zcorrelation[(int) (pos_min) -1];
	  					y_q = cmin;
	  					y_r = zcorrelation[(int) (pos_min) +1];
                        float denom = y_r+y_p-2*y_q;
						x_min = (denom>EPS*cmin)? (y_p-y_r)/(2*denom): 0.0f;
	  				
						subpixeldisparity[l] = disparity[l] + (pos_min+ x_min)/Z + T;
					}
				} 
		}
	}
	
	/*free memory*/
	free(zinput1);
	free(zinput2);
	free(correlation);
	free(zcorrelation);
}
