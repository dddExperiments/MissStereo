#define _USE_MATH_DEFINES // For Windows
#include "fft.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>

#define SIZE_PRIME 460
#define MAX_PRIME 3257


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr;


static int primes[SIZE_PRIME] = {2,3,5,7,11,13,17,19,23,29
,31,37,41,43,47,53,59,61,67,71
,73,79,83,89,97,101,103,107,109,113
,127,131,137,139,149,151,157,163,167,173
,179,181,191,193,197,199,211,223,227,229
,233,239,241,251,257,263,269,271,277,281
,283,293,307,311,313,317,331,337,347,349
,353,359,367,373,379,383,389,397,401,409
,419,421,431,433,439,443,449,457,461,463
,467,479,487,491,499,503,509,521,523,541
,547,557,563,569,571,577,587,593,599,601
,607,613,617,619,631,641,643,647,653,659
,661,673,677,683,691,701,709,719,727,733
,739,743,751,757,761,769,773,787,797,809
,811,821,823,827,829,839,853,857,859,863
,877,881,883,887,907,911,919,929,937,941
,947,953,967,971,977,983,991,997,1009,1013
,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069
,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151
,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223
,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291
,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373
,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451
,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511
,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583
,1597,1601,1607,1609,1613,1619,1621,1627,1637,1657
,1663,1667,1669,1693,1697,1699,1709,1721,1723,1733
,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811
,1823,1831,1847,1861,1867,1871,1873,1877,1879,1889
,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987
,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053
,2063,2069,2081,2083,2087,2089,2099,2111,2113,2129
,2131,2137,2141,2143,2153,2161,2179,2203,2207,2213
,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287
,2293,2297,2309,2311,2333,2339,2341,2347,2351,2357
,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423
,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531
,2539,2543,2549,2551,2557,2579,2591,2593,2609,2617
,2621,2633,2647,2657,2659,2663,2671,2677,2683,2687
,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741
,2749,2753,2767,2777,2789,2791,2797,2801,2803,2819
,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903
,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999
,3001,3011,3019,3023,3037,3041,3049,3061,3067,3079
,3083,3089,3109,3119,3121,3137,3163,3167,3169,3181
,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257};



void  md_fsig_copy(float *u,float *v,int size) { int i; for( i=0; i<size ;i++)  v[i]=u[i]; }


void md_fsig_clear(float *u,float value,int size) { int i; for( i=0; i<size; i++) u[i]=value;  }




/* decompose n into prime factors and returns the number of terms */
/* tab should be large enough to contain all factors (32 seems enough) */
/* Note: returns 0 if n==1, and dosn't work if n > MAX_PRIME^2 */

static int decompose(int * tab, int n)
{
  int count,i,p;

  if (n==1) return 0;

  /* search factors */
  for (count=i=0;i<SIZE_PRIME;i++)
    if ((n % primes[i])==0) {
      p=primes[i];
      do {
	tab[count]=p;
	count++; 
	n=n/p; 
      } while ((n%p)==0);
    }
  if (n!=1) { /* n is prime */
    tab[count]=n;
    count++;
  }
	
  return count;
}

void old_fft1d(float *Xr,float* Xi,float *Yr,float * Yi,int inverse, int size)
{
 
 
  
  int n,mmax,m,j,istep,i,isign;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  float *data;
  

  if (!inverse) isign=1; else isign=-1;
 
  /*float *data=new float[2*size];*/
  data = (float *)malloc(2*size*sizeof(float));
  
  for (i=0;i<size;i++){
  
  	data[2*i]=Xr[i];
    	data[2*i+1]=(Xi?Xi[i]:0.0f);
    }
  
  /*--- compute FFT of "data" array ---*/
  n=size << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j>i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    
    m=n >> 1;
    while (m >=2 && j>m) {
      j=j-m;
      m >>=1;
    }
    j=j+m;
  }
  
  mmax=2;
  while (n>mmax) {
    istep=2*mmax;
    theta=2*M_PI/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr=-2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m-1;i<=n-1;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i]=data[i]+tempr;
	data[i+1]=data[i+1]+tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }


   /* fill output */
    if (Yr) {
      if (inverse) 
	for (i=0;i<size;i++)
	  Yr[i]=data[2*i]/(float)size;
      else
	for (i=0;i<size;i++)
	  Yr[i]=data[2*i];
    }
  
  
  
    if (Yi) {
    
      if (inverse) 
	for (i=0;i<size;i++)
	  Yi[i]=data[2*i+1]/(float)size;
      else 
	for (i=0;i<size;i++)
	  Yi[i]=data[2*i+1];
    }
    
    free(data);
    
}



void fft1d(float *Xr,float *Xi,float *Yr,float *Yi, int inverse, int n)
{

  int tab[32],t,i,j,k,e,l,p,m,idxd,idxs,nsmp,mp;
  float *tra,*s,*d,*tmpf,wljx,wljy,*mc,*ms,isign;
  
  
  /*Decompose n*/
  t=decompose(tab,n);
  
  
  if (n>1 && tab[t-1]!=2) 
  {
  
    tra=(float *)malloc(4*n*sizeof(float));
    mc=(float *)malloc(n*sizeof(float));
    ms=(float *)malloc(n*sizeof(float));
   
    if (!inverse) isign=1; else isign=-1;

    for (i=0;i<n;i++) {
      mc[i]=cos(2*M_PI*((float)i)/((float)n));
      ms[i]=isign*sin(2*M_PI*((float)i)/((float)n));
    }

    s=tra;
    d=tra+2*n;

    
    /* fill source */
    if (Xi ==NULL)
      for(i=0;i<n;i++) {
	s[2*i]=Xr[i];
	s[2*i+1]=0;
      }
    else 
      for(i=0;i<n;i++) {
	s[2*i]=Xr[i];
	s[2*i+1]=Xi[i];
      }




    m=1;
    for (e=0;e<t;e++) {
      p=tab[e];
      nsmp=n/m/p;
      mp=m*p;
      for (k=0;k<2*n;k++) d[k]=0;
      for(j=0;j<p;j++) {
	for(l=0;l<mp;l++) {
	  wljx=mc[((l*j)%(mp))*(nsmp)];
	  wljy=ms[((l*j)%(mp))*(nsmp)];
	  for(i=0;i<nsmp;i++) {
	    idxd=i+l*nsmp;
	    idxs=i+(j+(l%m)*p)*(nsmp);
	    d[2*idxd]  +=  s[2*idxs]*wljx- s[2*idxs+1]*wljy;
	    d[2*idxd+1]+=  s[2*idxs]*wljy+ s[2*idxs+1]*wljx;
	  }
	}
      }
      tmpf=s;
      s=d;
      d=tmpf;
      m*=p;
    }
  
  

    /* fill output */
    if (Yr) {
      if (inverse) 
	for (i=0;i<n;i++)
	  Yr[i]=s[2*i]/(float)n;
      else
	for (i=0;i<n;i++)
	  Yr[i]=s[2*i];
    }
  
  
  
    if (Yi) {
    
      if (inverse) 
	for (i=0;i<n;i++)
	  Yi[i]=s[2*i+1]/(float)n;
      else 
	for (i=0;i<n;i++)
	  Yi[i]=s[2*i+1];
    }
    
    free(ms);
    free(mc);
    free(tra);
  
  } else 
  {
  
   old_fft1d(Xr,Xi,Yr,Yi,inverse,n);
  
  
  }
  
    
}







void fft1dzoom(float *in, float *out, int simetry, int zoom, int n)
{


  float *in_re,*in_im,*out_re,*out_im;
  float *in_part , *out_part;
  int nxz,hx,i,j;
   
  
  if (simetry) n = 2*n;
  
  in_re = (float *) malloc(n * sizeof(float));
  in_im = (float *) malloc(n * sizeof(float));
  in_part = (float *) malloc(n * sizeof(float));
  

  if (simetry) /*mirror simmetry to signal*/
  {
  
	  for(i=0;i<n;i++){

	    if (i<n/2) in_part[i]=in[n/2 -i-1];
	    else  in_part[i]=in[i- n/2];

	  }
	  
  } else 
  {
  
  	for(i = 0; i < n; i++) in_part[i] = in[i];
  }
  
  
  
  fft1d(in_part,NULL,in_re,in_im,0,n);



  nxz = n*zoom;

  out_re = (float *) malloc(nxz * sizeof(float));
  out_im = (float *) malloc(nxz * sizeof(float));
  out_part = (float *) malloc(nxz * sizeof(float));


  for(i=0; i < nxz ; i++)
  {
  	out_re[i] = 0.0;
  	out_im[i] = 0.0;
  }

  

  hx = n/2;
  for(i=0; i < hx + 1 ; i++){
    out_re[i] = in_re[i];
    out_im[i] = in_im[i];
  }


  for(; i < hx + 1 + (zoom-1) * n ; i++){
    out_re[i] = 0.0;
    out_im[i] = 0.0;
  }


  j= hx + 1;	
  for(; j < n; j++,i++){
    out_re[i] = in_re[j];
    out_im[i] = in_im[j];
  }
  
  
 fft1d(out_re,out_im,out_part,NULL,1,nxz);


  if (simetry) /*mirror simmetry to signal*/
  {
  
  	  for(i=0;i<nxz/2;i++)
	    out[i]= out_part[nxz/2 + i];

  
  } else 
  {
  
  	for(i = 0; i < nxz; i++) out[i] = out_part[i];
  }
  


 free(in_re);
 free(in_im);
 free(out_re);
 free(out_im);
 free(in_part);
 free(out_part);

}





/*in_re, out_re and out_img has memory*/
void fft2d(float *in_re, float *in_im,float *out_re, float *out_im, int i_flag, int width, int height)
{
  int      x,y,nx,ny, size;
  float  *f1_re,*f1_im,*f2_re,*f2_im;
  /*float    *dre_flag,*dim_flag;*/
  
  
  if ((!out_re) && (!out_im)) 
  {
  	return;
  }
  
    
  ny = height;
  nx = width;
  size = width * height;
  
  
  md_fsig_copy(in_re,out_re,size);
    
  if (in_im)
     md_fsig_copy(in_im,out_im,size);
  else 
     md_fsig_clear(out_im,0.0,size);

        
  /***** 1D FFT on lines *****/

f1_re = (float *) malloc(nx * sizeof(float));
f1_im = (float *) malloc(nx * sizeof(float));
f2_re = (float *) malloc(nx * sizeof(float));
f2_im = (float *) malloc(nx * sizeof(float));

    
  for (y=0;y<ny;y++) {
    
    for (x=0;x<nx;x++) {
      f1_re[x] = out_re[nx*y+x];
      f1_im[x] = out_im[nx*y+x];
    }
    
    fft1d(f1_re,f1_im,f2_re,f2_im,i_flag,nx);
    
    for (x=0;x<nx;x++) {
      out_re[nx*y+x] = f2_re[x];
      out_im[nx*y+x] = f2_im[x]; 
    }
  
  }
  
  free(f1_re);
  free(f1_im);
  free(f2_re);
  free(f2_im);

  /***** 1D FFT on columns *****/
  f1_re = (float *) malloc(ny * sizeof(float));
  f1_im = (float *) malloc(ny * sizeof(float));
  f2_re = (float *) malloc(ny * sizeof(float));
  f2_im = (float *) malloc(ny * sizeof(float));
 
  for (x=0;x<nx;x++) {
  
    for (y=0;y<ny;y++) {
      f1_re[y] = out_re[nx*y+x];
      f1_im[y] = out_im[nx*y+x];
    }
  
    fft1d(f1_re,f1_im,f2_re,f2_im,i_flag,ny);
  
    for (y=0;y<ny;y++) {
      out_re[nx*y+x] = f2_re[y];
      out_im[nx*y+x] = f2_im[y]; 
    }
  
  }	

  free(f1_re);
  free(f1_im);
  free(f2_re);
  free(f2_im);
  
}




void fftzoom(float *in,float *out,float z, int width, int height)
{
  float *in_re,*in_im,*out_re,*out_im,*aux;
  int nx,ny,nxz,nyz,x,y,adr,adrz,size,zsize;
  float z2,factor,zoom;


  size = width * height;
  nx = width;
  ny = height;

  zoom = z;
  if (zoom<=0.) { printf("Zoom factor must be positive.\n"); return; }
  
  
  /* DIRECT FFT */
  in_re = (float *) malloc(size * sizeof(float));
  in_im = (float *) malloc(size * sizeof(float));
    
  fft2d(in,NULL,in_re,in_im,0,width,height);


  /*for(int ii=0; ii < size; ii++) printf("%d: %f \t %f\n",ii,in_re[ii],in_im[ii]);*/

  
  /* COMPUTE DIMENSIONS */ 
  nxz = (int)floor((double)nx*(double)zoom+0.5);
  nyz = (int)floor((double)ny*(double)zoom+0.5);
  
  z2 = zoom*zoom;
  
  /* ALLOCATE AND INITIALIZE IMAGES */
  zsize = nxz * nyz;
  
  out_re = (float *) malloc(zsize * sizeof(float));
  out_im = (float *) malloc(zsize * sizeof(float));

  md_fsig_clear(out_re,0.0,zsize);
  md_fsig_clear(out_im,0.0,zsize);
  
  /* DISPATCH FOURIER COEFFICIENTS */
  if (zoom<1.) {
    
    /* UNZOOM */

    for (x=-nxz/2;x<=nxz/2;x++)
      for (y=-nyz/2;y<=nyz/2;y++) {
        adr  = (x+nx )%nx +((y+ny )%ny )*nx ;
        adrz = (x+nxz)%nxz+((y+nyz)%nyz)*nxz;
        out_re[adrz] += in_re[adr] * z2;
        out_im[adrz] += in_im[adr] * z2;
      }
  
  } else {
    /* ZOOM */
    
    for (x=-nx/2;x<=nx/2;x++)
      for (y=-ny/2;y<=ny/2;y++) {
        adr  = (x+nx )%nx +((y+ny )%ny )*nx ;
        adrz = (x+nxz)%nxz+((y+nyz)%nyz)*nxz;
        factor = z2 
	  * ((2*x==nx || 2*x==-nx)?0.5f:1.0f) 
	  * ((2*y==ny || 2*y==-ny)?0.5f:1.0f);
        out_re[adrz] = in_re[adr] * factor;
        out_im[adrz] = in_im[adr] * factor;
      }
  }
  
  /*float *aux = new float[zsize];*/
  aux=(float *)malloc(zsize*sizeof(float));

  md_fsig_clear(aux,0.0,zsize);
  
  /* INVERSE FFT */
  fft2d(out_re,out_im,out,aux,1,nxz,nyz);
  
  free(aux);
  free(in_re);
  free(in_im);
  free(out_re);
  free(out_im);
  
}
