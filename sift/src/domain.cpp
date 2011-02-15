#include "domain.h"
#include "numerics.h"

#define DEBUG 0

void compute_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, float **H)
{

/////////////////////////////////////////////////////// Normalize points to have baricenter [0][0] and standard deviation sqrt(2.0f)


	//if (DEBUG) printf("Compute homography with %d matches\n", n);


////////////////// Compute Baricenter
	float lx = 0.0,  ly = 0.0;
	float rx = 0.0,  ry = 0.0;

	for (int i=0; i< n; i++) {

		lx += x0[i]; 
		ly += y0[i]; 

		rx += x1[i]; 
		ry += y1[i]; 
	}

	lx /= (float) n; 
	ly /= (float) n; 
	rx /= (float) n; 
	ry /= (float) n; 


/////////////// Normalize points without modifying original vectors

	float *px0 = new float[n];
	float *py0 = new float[n];

	float *px1 = new float[n];
	float *py1 = new float[n];

	float spl= 0.0f, spr = 0.0f;
	for(int i=0; i < n; i++)
	{

		px0[i] = x0[i] - lx;
		py0[i] = y0[i] - ly;

		px1[i] = x1[i] - rx;
		py1[i] = y1[i] - ry;

		spl += sqrtf(px0[i] * px0[i] + py0[i]*py0[i]);
		spr += sqrtf(px1[i] * px1[i] + py1[i]*py1[i]);

	}


	spl = sqrtf(2.0f) / spl;
	spr = sqrtf(2.0f) / spr;

	for (int i=0; i< n; i++) {

		px0[i] *= spl;
		py0[i] *= spl;

		px1[i] *= spr;
		py1[i] *= spr;

	}

	
/////////////////////////////////////////////////////// Minimization problem || Ah ||


	float **Tpl = allocate_float_matrix(3,3);
	float **Timinv = allocate_float_matrix(3,3);

	// similarity transformation of the plane 
	Tpl[0][0] = spl; Tpl[0][1] = 0.0; Tpl[0][2] = -spl*lx;
	Tpl[1][0] = 0.0; Tpl[1][1] = spl; Tpl[1][2] = -spl*ly;
	Tpl[2][0] = 0.0; Tpl[2][1] = 0.0; Tpl[2][2] = 1.0;
	
	// inverse similarity transformation of the image
	Timinv[0][0] = 1.0/spr; Timinv[0][1] =   0.0  ; Timinv[0][2] = rx;
	Timinv[1][0] =   0.0  ; Timinv[1][1] = 1.0/spr; Timinv[1][2] = ry;
	Timinv[2][0] =   0.0  ; Timinv[2][1] =   0.0  ; Timinv[2][2] = 1.0;


///////////////// System matrix

	float **A = allocate_float_matrix(2*n,9);
	for(int i=0, eq=0; i < n; i++, eq++) {

		float	xpl = px0[i], ypl = py0[i],
			xim = px1[i], yim = py1[i];

		A[eq][0] = A[eq][1] = A[eq][2] = 0.0;
		A[eq][3] = -xpl;
		A[eq][4] = -ypl;
		A[eq][5] = -1.0;
		A[eq][6] =  yim * xpl;
		A[eq][7] =  yim * ypl;
		A[eq][8] =  yim;

		eq++;

		A[eq][0] =  xpl;
		A[eq][1] =  ypl;
		A[eq][2] =  1.0;
		A[eq][3] = A[eq][4] = A[eq][5] = 0.0;
		A[eq][6] = -xim * xpl;
		A[eq][7] = -xim * ypl;
		A[eq][8] = -xim;
	}


///////////////// SVD
	float **U = allocate_float_matrix(2*n,9);
	float **V = allocate_float_matrix(9,9);
	float *W = new float[9];

	compute_svd(A,U,V,W,2*n,9);

	
	// Find the index of the least singular value
	int imin = 0;
	for (int i=1; i < 9; i++) 
		if ( W[i] < W[imin] ) imin = i;


////////////////// Denormalize H = Timinv * V.col(imin)* Tpl;
	float **matrix = allocate_float_matrix(3,3);
	float **result = allocate_float_matrix(3,3);


	int k=0;
	for(int i=0; i < 3; i++)
		for(int j=0; j < 3; j++, k++)
			matrix[i][j] = V[k][imin];

	
	product_square_float_matrixes(result,Timinv,matrix,3);

	product_square_float_matrixes(H,result,Tpl,3);


	desallocate_float_matrix(U ,2*n,9);
	desallocate_float_matrix(V,9,9);
	delete[] W;
	
	desallocate_float_matrix(Tpl,3,3);
	desallocate_float_matrix(Timinv,3,3);
	desallocate_float_matrix(matrix,3,3);
	desallocate_float_matrix(result,3,3);
	desallocate_float_matrix(A,2*n,9);

	delete[] px0;
	delete[] py0;
	delete[] px1;
	delete[] py1;


}





/// Compute homography using svd + Ransac
void compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, float **H)
{

	// Initialize seed		
	srand( (long int) time (NULL));	
	float tolerance2 = tolerance * tolerance;


	float **Haux = allocate_float_matrix(3,3);

	int *accorded = new int[n];	
	int *paccorded = new int[n];	
	int naccorded = 0;
	int pnaccorded = 0;

	for(int iter = 0; iter < niter; iter++)
	{


		// Choose 4 indexos from 1..n without repeated values
		int indexos[4];
		int acceptable = 0;
		while (!acceptable)
		{
			acceptable = 1;
			for(int i=0; i < 4; i++) indexos[i] = (int)  floor(rand()/(double)RAND_MAX * (double) n);

			// Check if indexos are repeated
			for(int i=0; i < 4; i++)
				for(int j=i+1; j < 4; j++)
					if (indexos[i] == indexos[j]) acceptable = 0; 
		}


		


		// Store selected matches 
		float px0[4] , py0[4], px1[4] , py1[4];
		for(int i=0; i < 4; i++)
		{
			px0[i] = x0[indexos[i]];
			py0[i] = y0[indexos[i]];

			px1[i] = x1[indexos[i]];
			py1[i] = y1[indexos[i]];
		}


		// Compute planar homography
		compute_planar_homography_n_points(px0, py0, px1, py1, 4, Haux);


		// Which matches are according to this transformation
		pnaccorded = 0;
		for(int i=0; i < n; i++)
		{

			float vec[3];
			vec[0] = x0[i];
			vec[1] = y0[i];
			vec[2] = 1.0f;
		
			float res[3];
			float_vector_matrix_product(Haux, vec ,res , 3);

			if (res[2] != 0.0f) {

				res[0] /= res[2]; res[1] /= res[2];

				float dif = (res[0] - x1[i]) * (res[0] - x1[i]) + (res[1] - y1[i]) * (res[1] - y1[i]);
			
				if (dif < tolerance2) {  paccorded[pnaccorded] = i; pnaccorded++; }
			
			}
			
		}


		if (DEBUG) printf("%d: %d %d %d %d --> %d \n", iter, indexos[0],indexos[1],indexos[2],indexos[3], pnaccorded);

		
		if (pnaccorded > naccorded)
		{

			naccorded = pnaccorded;
			for(int i=0; i < naccorded; i++) accorded[i] = paccorded[i];
			for(int i=0; i < 3; i++) for(int j=0; j < 3; j++) H[i][j] = Haux[i][j];

		}

	}


	desallocate_float_matrix(Haux,3,3);

}




/* logarithm (base 10) of binomial coefficient */
float logcombi(int k,int n)
{
  double r;
  int i;
 
  if (k>=n || k<=0) return(0.);
  if (n-k<k) k=n-k;
  r = 0.;
  for (i=1;i<=k;i++) 
    r += log10((double)(n-i+1))-log10((double)i);

  return((float)r);
}

/* tabulate logcombi(.,n) */
float * makelogcombi_n(int n)
{
  float *l;
  int k;

  l = (float *)malloc((n+1)*sizeof(float));
  for (k=0;k<=n;k++) l[k]=logcombi(k,n);

  return(l);
}

/* tabulate logcombi(k,.) */
float *makelogcombi_k(int k,int nmax)
{
  float *l;
  int n;

  l = (float *)malloc((nmax+1)*sizeof(float));
  for (n=0;n<=nmax;n++) l[n]=logcombi(k,n);

  return(l);
}





void compute_moisan_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float &epsilon, float **H, int &counter, int recursivity)
{



	// Initialize seed		
	srand( (long int) time (NULL) );	
	

	float **Haux = allocate_float_matrix(3,3);

	float *dist = new float[n];
	float *dindexos = new float[n];

	int nselected = 0;
	float *i0selected = new float[n];
	float *j0selected = new float[n];

	float *i1selected = new float[n];
	float *j1selected = new float[n];


  	/* tabulate logcombi */
	float loge0 = (float)log10((double)(n-4));
  	float *logcn = makelogcombi_n(n);
	float *logc4 = makelogcombi_k(4,n); 


	// Normalize points using minimum and maximum coordinates 
	float xmin, xmax, ymin, ymax;
	xmin = xmax = x1[0];
	ymin = ymax = x1[0];

	for(int i=0; i < n; i++)
	{
		if (x1[i] < xmin) xmin = x1[i];
		if (y1[i] < ymin) ymin = y1[i];

		if (x1[i] > xmax) xmax = x1[i];
		if (y1[i] > ymax) ymax = y1[i];


	}

	float mepsilon = 10000000.0f;
	for(int iter = 0; iter < niter; iter++)
	{

		// Initializing distances and indexos for the whole vector
		for(int i=0; i < n ; i++)
		{

			dist[i] = 0.0;
			dindexos[i] = (float) i;
		}



		// Choose 4 indexos from 1..n without repeated values
		int indexos[4];
		int acceptable = 0;
		while (!acceptable)
		{
			acceptable = 1;
			for(int i=0; i < 4; i++) indexos[i] = (int)  floor(rand()/(double)RAND_MAX * (double) n);

			// Check if indexos are repeated
			for(int i=0; i < 4; i++)
				for(int j=i+1; j < 4; j++)
					if (indexos[i] == indexos[j]) acceptable = 0; 
		}


		


		// Store selected matches 
		float px0[4] , py0[4], px1[4] , py1[4];
		for(int i=0; i < 4; i++)
		{
			px0[i] = x0[indexos[i]];
			py0[i] = y0[indexos[i]];

			px1[i] = x1[indexos[i]];
			py1[i] = y1[indexos[i]];
		}


		// Compute planar homography
		compute_planar_homography_n_points(px0, py0, px1, py1, 4, Haux);



		for(int i=0; i < n; i++)
		{

			float vec[3];
			vec[0] = x0[i];
			vec[1] = y0[i];
			vec[2] = 1.0f;
		
			float res[3];
			float_vector_matrix_product(Haux, vec ,res , 3);

			if (res[2] != 0.0f) {

				res[0] /= res[2]; res[1] /= res[2];

				float dif = (res[0] - x1[i]) * (res[0] - x1[i]) + (res[1] - y1[i]) * (res[1] - y1[i]);
			
				dist[i] = dif;

				//if (dif < tolerance2) {  paccorded[pnaccorded] = i; pnaccorded++; }
			
			} else dist[i] = 100000000.0f;

			
		}


		// Order distances
		quick_sort(dist,dindexos,n);	
	

		// Look for most meaningful subset
		for(int i=4 ; i < n; i++)
		{

			float rigidity = dist[i];	
			float logalpha = 0.5f*(float)log10((double) rigidity);
					
			float nfa = (float)  loge0  +  (float)(i-3) * logalpha + logcn[i+1] + logc4[i+1] ;

			if (nfa < mepsilon) 
			{

				mepsilon = nfa;
				for(int ki=0; ki < 3; ki++) for(int kj=0; kj < 3; kj++) H[ki][kj] = Haux[ki][kj];

				nselected = i;
	
				for(int ki=0; ki < nselected; ki++) {
						i0selected[ki] = x0[(int) dindexos[ki]];
						j0selected[ki] = y0[(int) dindexos[ki]];
						i1selected[ki] = x1[(int) dindexos[ki]];
						j1selected[ki] = y1[(int) dindexos[ki]];
				}
			
							
			}

		}


		
	}


	epsilon = mepsilon;
	counter++;

	if (counter < recursivity) 								
		compute_moisan_planar_homography_n_points(i0selected,j0selected,i1selected,j1selected,nselected,niter,epsilon,H,counter,recursivity);



	desallocate_float_matrix(Haux,3,3);
	delete[] dist;
	delete[] dindexos;

	delete[] i0selected;
	delete[] j0selected;
	delete[] i1selected;
	delete[] j1selected;

}




void compute_planar_homography_bounding_box(int width, int height, float **H, float *x0, float *y0,  int *nwidth, int *nheight)
{


	int xmin, xmax, ymin , ymax;
	ymin = xmin = 10000000;
	ymax = xmax = -10000000;

	for(int i = 0; i <= 1; i++)
		for(int j = 0; j <= 1; j++)
	{

		float vec[3], res[3];

		vec[0] = (float) (i*(width)); 
		vec[1] = (float) (j*(height)); 
		vec[2] = 1.0f; 

		float_vector_matrix_product(H, vec ,res , 3);


		if (res[2] != 0.0f)
		{
			res[0] /= res[2]; res[1] /= res[2];
			
			xmin = MIN(xmin , res[0]);	xmax = MAX(xmax , res[0]);
			ymin = MIN(ymin , res[1]);	ymax = MAX(ymax , res[1]);

		}

	}


	*nwidth = (int) ceilf(xmax - xmin);
	*nheight = (int) ceilf(ymax - ymin);
	*x0 = xmin;
	*y0 = ymin;
}



void apply_planar_homography(float *input, int width, int height, float **H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight)
{


	if (DEBUG) printf("width: %d height: %d ------>  nwidth: %d  nheight: %d \n", width, height, nwidth, nheight);	

	/// We compute inverse transformation

	if (DEBUG) 
	{

		printf("Matrix: \n");
		print_float_matrix(H,3,3);
	}
	
	float **V = allocate_float_matrix(3,3);

 	luinv(H, V, 3);

	if (DEBUG) 
	{

		printf("Inverse Matrix: \n");
		print_float_matrix(V,3,3);
	}



	float *coeffs;
	float *ref;

	float  cx[12],cy[12],ak[13];

	if (order!=0 && order!=1 && order!=-3 && 
      order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
    	{	
		printf("unrecognized interpolation order.\n");
		exit(-1);
	}

        if (order>=3) {

		coeffs = new float[width*height];
    		finvspline(input,order,coeffs,width,height);
    		ref = coeffs;
    		if (order>3) init_splinen(ak,order);

	} else 
	{
    		coeffs = NULL;
    		ref = input;
  	}

	int xi,yi;
	float xp,yp;
	float res;
	int n1,n2;
	float p=-0.5;

	/// For each point in new image we compute its anti image and interpolate the new value
	for(int i=0; i < nwidth; i++)
		for(int j=0; j < nheight; j++)
	{


		float vec[3];
		vec[0] = (float) i + x0;
		vec[1] = (float) j + y0;
		vec[2] = 1.0f;


		float vres[3];
		float_vector_matrix_product(V, vec ,vres , 3);

		if (vres[2] != 0.0f) 
		{

			vres[0] /= vres[2]; vres[1] /= vres[2];
				
			

			xp =  (float) vres[0];
			yp =  (float) vres[1];

			if (order == 0) { 
	
				xi = (int)floor((double)xp); 
				yi = (int)floor((double)yp);
		
				if (xi<0 || xi>=width || yi<0 || yi>=height)
		 			 res = bg; 
				else res = input[yi*width+xi];
	
     			 } else { 
	
		
				if (xp<0. || xp>=(float)width || yp<0. || yp>=(float)height) res=bg; 
				else {
					//xp -= 0.5; yp -= 0.5;
					int xi = (int)floor((double)xp); 
					int yi = (int)floor((double)yp);
					float ux = xp-(float)xi;
	  				float uy = yp-(float)yi;

					switch (order) 
	   				{
					    	case 1: /* first order interpolation (bilinear) */
	      						n2 = 1;
							cx[0]=ux;	cx[1]=1.f-ux;
							cy[0]=uy; cy[1]=1.f-uy;
							break;
						
						case -3: /* third order interpolation (bicubic Keys' function) */
							n2 = 2;
							keys(cx,ux,p);
							keys(cy,uy,p);
							break;

						case 3: /* spline of order 3 */
							n2 = 2;
							spline3(cx,ux);
							spline3(cy,uy);
							break;

						default: /* spline of order >3 */
							n2 = (1+order)/2;
							splinen(cx,ux,ak,order);
							splinen(cy,uy,ak,order);
							break;
					}
	  
	  				res = 0.; n1 = 1-n2;
	 				if (xi+n1>=0 && xi+n2<width && yi+n1>=0 && yi+n2<height) {
	    				
						int adr = yi*width+xi; 
	    					for (int dy=n1;dy<=n2;dy++) 
	    					  for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*ref[adr+width*dy+dx];
	  				} else 
	
					   	for (int dy=n1;dy<=n2;dy++)
						      for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,bg,width,height);
					
      				}

			}	

		

			out[j*nwidth+i] = res;
   			

		}
	
	}


}






void apply_zoom(float *input, float *out, float zoom, int order, int width, int height)
{

	int nwidth = (int)( zoom * (float) width);
	int nheight = (int)( zoom * (float) height);

	float *coeffs;
	float *ref;

	float  cx[12],cy[12],ak[13];

	if (order!=0 && order!=1 && order!=-3 && 
      order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
    	{	
		printf("unrecognized interpolation order.\n");
		exit(-1);
	}

        if (order>=3) {

		coeffs = new float[width*height];
    		finvspline(input,order,coeffs,width,height);
    		ref = coeffs;
    		if (order>3) init_splinen(ak,order);

	} else 
	{
    		coeffs = NULL;
    		ref = input;
  	}

	int xi,yi;
	float xp,yp;
	float res;
	int n1,n2;
	float bg = 0.0f;
	float p=-0.5;
	for(int i=0; i < nwidth; i++)
		for(int j=0; j < nheight; j++)
		{

			xp =  (float) i / zoom;
			yp =  (float) j / zoom;

			if (order == 0) { 
	
				xi = (int)floor((double)xp); 
				yi = (int)floor((double)yp);
		
				if (xi<0 || xi>=width || yi<0 || yi>=height)
		 			 res = bg; 
				else res = input[yi*width+xi];
	
     			 } else { 
	
		
				if (xp<0. || xp>=(float)width || yp<0. || yp>=(float)height) res=bg; 
				else {
					xp -= 0.5; yp -= 0.5;
					int xi = (int)floor((double)xp); 
					int yi = (int)floor((double)yp);
					float ux = xp-(float)xi;
	  				float uy = yp-(float)yi;

					switch (order) 
	   				{
					    	case 1: /* first order interpolation (bilinear) */
	      						n2 = 1;
							cx[0]=ux;	cx[1]=1.f-ux;
							cy[0]=uy; cy[1]=1.f-uy;
							break;
						
						case -3: /* third order interpolation (bicubic Keys' function) */
							n2 = 2;
							keys(cx,ux,p);
							keys(cy,uy,p);
							break;

						case 3: /* spline of order 3 */
							n2 = 2;
							spline3(cx,ux);
							spline3(cy,uy);
							break;

						default: /* spline of order >3 */
							n2 = (1+order)/2;
							splinen(cx,ux,ak,order);
							splinen(cy,uy,ak,order);
							break;
					}
	  
	  				res = 0.; n1 = 1-n2;
	 				if (xi+n1>=0 && xi+n2<width && yi+n1>=0 && yi+n2<height) {
	    				
						int adr = yi*width+xi; 
	    					for (int dy=n1;dy<=n2;dy++) 
	    					  for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*ref[adr+width*dy+dx];
	  				} else 
	
					   	for (int dy=n1;dy<=n2;dy++)
						      for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,bg,width,height);
					
      			}

		}	

		out[j*nwidth+i] = res;
    	
	}
    if(ref != input)
        delete [] ref;
}



void apply_general_transformation(float *input, float *transformx, float* transformy, float *out, float bg, int order, int width, int height, int nwidth,int nheight)
{

	float *coeffs;
	float *ref;

	float  cx[12],cy[12],ak[13];

	if (order!=0 && order!=1 && order!=-3 && 
      order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
    	{	
		printf("unrecognized interpolation order.\n");
		exit(-1);
	}

        if (order>=3) {

		coeffs = new float[width*height];
    		finvspline(input,order,coeffs,width,height);
    		ref = coeffs;
    		if (order>3) init_splinen(ak,order);

	} else 
	{
    		coeffs = NULL;
    		ref = input;
  	}

	int xi,yi;
	float xp,yp;
	float res;
	int n1,n2;
	float p=-0.5;
	for(int i=0; i < nwidth; i++)
		for(int j=0; j < nheight; j++)
		{

			xp =  (float) transformx[j*nwidth+i];
			yp =  (float) transformy[j*nwidth+i];

			if (order == 0) { 
	
				xi = (int)floor((double)xp); 
				yi = (int)floor((double)yp);
		
				if (xi<0 || xi>=width || yi<0 || yi>=height)
		 			 res = bg; 
				else res = input[yi*width+xi];
	
			} else { 
	
		
				if (xp<0. || xp>(float)width || yp<0. || yp>(float)height) res=bg; 
				else {
					//xp -= 0.5; yp -= 0.5;
					
					int xi = (int)floor((double)xp); 
					int yi = (int)floor((double)yp);
					float ux = xp-(float)xi;
	  				float uy = yp-(float)yi;

					switch (order) 
	   				{
					    	case 1: /* first order interpolation (bilinear) */
	      						n2 = 1;
							cx[0]=ux;	cx[1]=1.f-ux;
							cy[0]=uy; cy[1]=1.f-uy;
							break;
						
						case -3: /* third order interpolation (bicubic Keys' function) */
							n2 = 2;
							keys(cx,ux,p);
							keys(cy,uy,p);
							break;

						case 3: /* spline of order 3 */
							n2 = 2;
							spline3(cx,ux);
							spline3(cy,uy);
							break;

						default: /* spline of order >3 */
							n2 = (1+order)/2;
							splinen(cx,ux,ak,order);
							splinen(cy,uy,ak,order);
							break;
					}
	  
	  				res = 0.; n1 = 1-n2;
	 				if (xi+n1>=0 && xi+n2<width && yi+n1>=0 && yi+n2<height) {
	    				
						int adr = yi*width+xi; 
	    					for (int dy=n1;dy<=n2;dy++) 
	    					  for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*ref[adr+width*dy+dx];
	  				} else 
	
					   	for (int dy=n1;dy<=n2;dy++)
						      for (int dx=n1;dx<=n2;dx++) 
							res += cy[n2-dy]*cx[n2-dx]*v(ref,xi+dx,yi+dy,bg,width,height);
					
      			}

		}	

		out[j*nwidth+i] = res;
    	
	}

}
