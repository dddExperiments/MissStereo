#include "stereoAC.h"
#include "libStereo/patch.h"
#include "libIO/nan.h"

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>

/// Probability quantification levels
static const int NB_PREC = 5;

/// Number of components to compare
static const int NB_COMPS = 12;

/// Interval of number of PCs tested for max criterion
static const int MIN_PC=5, MAX_PC=10;

/// Radius of window patch.
static int radius(int n) {
 	int r = (int)sqrtf((float)n);
	r = (r-1) / 2;
    assert((2*r+1)*(2*r+1) == n);
    return r;
}

/// Number of increasing maps from [1,npcused] to [1,nb_prec]
static float compute_inc_maps(int npcused, int nb_prec)
{
    float FC=0.0f;
    const int n = npcused+nb_prec-2;
    LWImage<float> C = alloc_image<float>(n,n); // Combinatorial coefficients

    for(int k=0; k<n; k++) {
        *C.pixel(0,k) = 1.0f;
        *C.pixel(k,k) = 1.0f;
    }

    for(int k=2; k<n; k++)
        for(int j=1; j<k; j++)
            *C.pixel(j,k) = *C.pixel(j-1,k-1) + *C.pixel(j,k-1);

    for(int k=0; k<nb_prec; k++)
        FC += (k+1) * *C.pixel(nb_prec-1-k,n-k-1);

    free(C.data);
    return FC;
}

/// Functor for comparison of indices in an array.
class CompareAbs {
    const float* vec;
public:
    CompareAbs(const float* v): vec(v) {}
    bool operator()(int i, int j) { // |vec[i]| > |vec[j]|
        float v = vec[i];
        if(v < 0) v = -v;
        return (v>vec[j] && v>-vec[j]);
    }
};

/// Return in \a ivector the indices of the \a size highest abs components in
/// \a vec.
static void max_components(const float* vec, int n,
                           std::vector<int>& ivector, int size)
{
    ivector.clear();
    for(int i=0; i<n; i++)
        ivector.push_back(i);
    std::partial_sort(ivector.begin(), ivector.begin()+size, ivector.end(),
                      CompareAbs(vec));
}

// Count valid patches (no NaN) in image
static int valid_patches(const LWImage<float>& im, int win)
{
    int nValid=0;
	for(int y=win; y+win < im.h; y++)
        for(int x=win; x+win < im.w; x++) {
            ++ nValid;
			for(int r=-win; r<=win; r++)
                for(int s=-win; s<=win; s++)
                    if(! is_number(*im.pixel(x+r,y+s))) {
                        r=s=win;
                        -- nValid;
                        continue;
                    }
        }
    return nValid;
}

// Compute descriptor of patch
static bool descriptor(const float* igray, int width, int win,
                       const LWImage<float>& PCs,
                       float* desc, int stride=1)
{
    const int n = PCs.w;
    for(int i=0; i<n; i++)
        desc[i*stride] = 0.0f;
    const float* pc = PCs.data;
    for(int r=-win; r<=win; r++)
        for(int s=-win; s<=win; s++) {
            float v = igray[s*width+r];
            for(int i=0; i<n; i++)
                desc[i*stride] += v * *pc++;
            if(! is_number(v))
                return false;
        }
    return true;
}

/// Coef: line=component, column=patch. Return number of valid patches.
static int compute_coefficients(const LWImage<float>& im, LWImage<float>& Coef,
                                const LWImage<float>& PCs, int win)
{
	int i=0, nValid=0;
    for(int y=0; y < im.h; y++)
        for(int x=0; x < im.w; x++)
            if(win<=x && x+win<im.w && win<=y && y+win<im.h) {
                if(descriptor(im.pixel(x,y), im.w, win, PCs,
                              Coef.pixel(i++,0), Coef.w))
                    ++nValid;
            } else { // Invalid patch at borders
                for(int j=0; j<Coef.h; j++)
                    *Coef.pixel(i,j) = NaN;
                ++i;
            }
    return nValid;
}

static bool is_valid(const std::pair<float,int>& p)
{ return is_number(p.first); }
static bool compare(const std::pair<float,int>& p1,
                    const std::pair<float,int>& p2)
{ return (p1.first < p2.first); }

/// Order each line of \a Coef independently. Each line of \a Index encodes the
/// applied permutation. Invalid patches (NaN) are gathered at the end, and
/// their index is negative.
static void order_coefficients(LWImage<float>& Coef, LWImage<int>& Index)
{
    std::pair<float,int>* aux = new std::pair<float,int>[Coef.w];

	for(int k=0; k < Coef.h; k++) {
		for(int j=0; j < Coef.w; j++)
            aux[j] = std::make_pair(*Coef.pixel(j,k),j);

        // Push NaN to the end
        std::pair<float,int>* it = std::partition(aux,aux+Coef.w,is_valid);
        // Sort only valid coefficients
        std::sort(aux, it, compare);

        const int end=it-aux;
		for(int j=0; j < Coef.w; j++) {
            *Coef.pixel(j,k) = aux[j].first;
            *Index.pixel(aux[j].second,k) = (j<end)? j: -1;
        }
	}

    delete [] aux;
}

/// Probability that a uniform random variable in [0,end2) be as close to
/// ind1 as ind2 is. \a norm should be 1/end2, and is used to avoid a division.
static float proba(int ind1, int ind2, int end2, float norm)
{
    int dif = ind2-ind1;
    if(dif < 0)
        dif = -dif;
    if(end2-ind1 < dif)
        return ((end2-ind2)*norm);
    else if(ind1 < dif)
        return (ind2*norm);
    return (2*dif*norm);
}

/// Quantify and make increasing each line of \a tab, according to quantization
/// levels in \a probaIni.
static void quantize_increase(LWImage<float>& tab,
                              const std::vector<float>& probaIni)
{
    std::vector<float>::const_iterator last = probaIni.end()-1;
    float* p=tab.data;
    for(int j=0; j<tab.h; j++) {
        std::vector<float>::const_iterator it = probaIni.begin();
        for(int i=0; i<tab.w; i++, p++) {
            while(it!=last && *it<*p)
                ++it;
            *p = *it;
        }
    }
}

/// Compute NFA from individual probabilities \a pac,
/// increasing sequence criterion.
static void compute_NFA_inc(const LWImage<float>& pac, float* NFA)
{
    const float* p = pac.data;
    for(int d=0; d<pac.h; d++) {
        float nfa = 1.0f;
        for(int i=0; i<pac.w; i++)
            nfa *= *p++;
        *NFA++ = nfa;
    }
}

/// Compute NFA from individual probabilities \a pac,
/// max criterion.
static void compute_NFA_max(const LWImage<float>& pac, float* NFA)
{
    for(int d=0; d<pac.h; d++) {
        float minNFA = FLT_MAX;
        const float* p = pac.pixel(MIN_PC,d);
        for(int m=MIN_PC; m<=MAX_PC; m++) {
            float nfa = pow(*p++, m+1);
            if(nfa < minNFA)
                minNFA = nfa;
        }
        *NFA++ = minNFA;
    }
}

/// Put in \a result the meaningful shift yielding lowest patch distance.
static void setBest(const float* NFA, int nDisp, int minDisp, float eps,
                    const LWImage<float>& im1, const LWImage<float>& im2,
                    int i, int j, int win, float& result)
{
    float bestSSD = FLT_MAX;
    for(int d=0; d<nDisp; d++)
        if(*NFA++ < eps) {
            float ssd = cssd(im1,i,j, im2,i+d+minDisp,j, win);
            if(ssd < bestSSD) {
                bestSSD = ssd;
                result = static_cast<float>(d+minDisp);
            }
        }
}

/// a contrario stereo: find meaningful patch matches from left to right.
/// \a PCs is the matrix projecting a patch into the principal components.
/// Disparities should be searched in interval [minDisp,maxDisp].
/// Two criteria are possible: increasing sequence and max (L_infinity).
void stereoAC(LWImage<float> im1, LWImage<float> im2, LWImage<float> PCs,
              int minDisp, int maxDisp, float epsNFA,
              float* dispInc, float* dispMax)
{
    const int win = radius(PCs.h);
    const int nComps = std::min(NB_COMPS, PCs.h);
    assert(MAX_PC < nComps);

  	// Compute coefficients
	LWImage<float> Coef2 = alloc_image<float>(im2.w*im2.h,PCs.h);
	const int end1 = valid_patches(im1, win);
    const int end2 = compute_coefficients(im2, Coef2, PCs, win);
    assert(end2 == valid_patches(im2, win));
    const float norm=1.0f/end2;

	// Order coefficients
  	LWImage<int> Index2 = alloc_image<int>(Coef2.w, Coef2.h);
	order_coefficients(Coef2, Index2);

	// Probability quantization: 1/2^(*nb_prec-1) ... 1/16 1/8 1/4 1/2 1
    std::vector<float> probaIni;
    probaIni.push_back(1.0f);
	for(int k=NB_PREC-2; k>=0; k--)
        probaIni.insert(probaIni.begin(), probaIni.front()/2.0f);

    // Number of tests
	const int nDisp = maxDisp-minDisp+1;
	float nTests = (float)end1 * (float)nDisp;
    const float epsInc = epsNFA/(nTests*compute_inc_maps(nComps,NB_PREC));
    const float epsMax = epsNFA/(nTests*(MAX_PC-MIN_PC+1));

    // Buffers
    float* desc = new float[PCs.w]; // Descriptor of patch
    std::vector<int> ivector; // Indices of highest components
    LWImage<float> pac = alloc_image<float>(nComps,nDisp); // Elementary proba
 	float* NFA = new float[nDisp];

    if(dispInc)
        std::fill_n(dispInc, im1.w*im1.h, NaN); // Init output image
    if(dispMax)
        std::fill_n(dispMax, im1.w*im1.h, NaN); // Init output image

    for(int y=win; y+win<im1.h && y<im2.h; y++)
        for(int x=win; x+win<im1.w; x++) {
            const int iPixel = y*im1.w+x;
            // Find components with max abs value in first image and order them
            if(! descriptor(&im1.data[iPixel], im1.w, win, PCs, desc))
                continue; // Invalid window: skip
            max_components(desc, PCs.w, ivector, nComps);

            std::fill_n(pac.data, pac.w*pac.h, 1.0f);
            for(int m=0; m<nComps; m++) {
                // Find place in ordered coefficients of second image
                int iComp = ivector[m]; // Component index
                const float* beg = Coef2.pixel(0,iComp);
                int ind1 = std::lower_bound(beg,beg+Coef2.w,desc[iComp]) - beg;

                int x2=x+minDisp, pos=y*im2.w+x2; //Patch number
                for(int d=0; d<nDisp; d++, x2++, pos++)
                    if(0<=x2 && x2<im2.w) {
                        int ind2 = *Index2.pixel(pos,iComp); // Place in histo
                        if(ind2 >= 0) // Valid patch
                            *pac.pixel(m,d) = proba(ind1, ind2, end2, norm);
                    }
            }

            quantize_increase(pac, probaIni);
            if(dispInc) {
                compute_NFA_inc(pac, NFA);
                setBest(NFA, nDisp, minDisp, epsInc, im1, im2, x, y, win,
                        dispInc[iPixel]);
            }
            if(dispMax) {
                compute_NFA_max(pac, NFA);
                setBest(NFA, nDisp, minDisp, epsMax, im1, im2, x, y, win,
                        dispMax[iPixel]);
            }
        }

	delete [] NFA;
    delete [] desc;

 	free(Coef2.data);
 	free(Index2.data);
 	free(pac.data);
}
