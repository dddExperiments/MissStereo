//
// C++ Implementation: stereomatch
//
// Description: 
//
//
// Author:  <>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "orsa.h"
#include "libNumerics/numerics.h"

#define _USE_MATH_DEFINES // For Windows (M_PI)
#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

// Logarithm (base 10) of binomial coefficient
static float logcombi(size_t k, size_t n)
{
    if (k>=n || k<=0) return(0.);
    if (n-k<k) k=n-k;
    double r = 0.;
    for(size_t i=1; i <= k; i++) 
        r += log10((double)(n-i+1))-log10((double)i);

    return (float)r;
}

// Tabulate logcombi(.,n)
static float* makelogcombi_n(size_t n)
{
    float* l = new float[n+1];
    for(size_t k=0; k<=n; k++)
        l[k] = logcombi(k,n);
    return l;
}

// Tabulate logcombi(k,.)
static float* makelogcombi_k(size_t k, size_t nmax)
{
    float* l = new float[nmax+1];
    for(size_t n=0; n <= nmax; n++)
        l[n]=logcombi(k,n);
    return l;
}

// Get a (sorted) random 7-uple of 0..n-1
static void random_p7(size_t k[7], size_t n)
{
    for(size_t i=0; i < 7; i++) {
        size_t r = (rand()>>3)%(n-i), j;
        for(j=0; j<i && r>=k[j]; j++)
            r++;
        size_t j0 = j;
        for(j=i; j > j0; j--)
            k[j]=k[j-1];
        k[j0] = r;
    }
}

// Return in k[0:6] the next sample of 7 values in range 0:k[7]-1
static bool next_sample(size_t k[8])
{ // next ordered sample
    int i0=6; // Find rightmost index we can increment
    while(i0>=0 && k[i0]+1==k[i0+1])
        --i0;
    if(i0 < 0) // No more sample
        return false;
    k[i0]++;
    for(int i=i0+1; i < 7; i++)
        k[i] = k[i-1]+1;
    return true;
}

// Compute the cubic polynomial det(F1+xF2) of variable x
static void characPoly(const libNumerics::matrix<float>& F1,
                       const libNumerics::matrix<float>& F2,
                       float a[4])
{
    assert(F1.nrow()==3 && F1.ncol()==3 && F2.nrow()==3 && F2.ncol()==3);
    a[0] = a[1] = a[2] = a[3] = 0.0f;
    for(int i0=0; i0 < 3; i0++) { // Even permutations
        int i1 = (i0+1)%3;
        int i2 = (i1+1)%3;
        a[0] += F1(i0,0)*F1(i1,1)*F1(i2,2);
        a[1] += F2(i0,0)*F1(i1,1)*F1(i2,2)+
                F1(i0,0)*F2(i1,1)*F1(i2,2)+
                F1(i0,0)*F1(i1,1)*F2(i2,2);
        a[2] += F1(i0,0)*F2(i1,1)*F2(i2,2)+
                F2(i0,0)*F1(i1,1)*F2(i2,2)+
                F2(i0,0)*F2(i1,1)*F1(i2,2);
        a[3] += F2(i0,0)*F2(i1,1)*F2(i2,2);
    }
    for(int i0=0; i0 < 3; i0++) { // Odd permutations
        int i1 = (i0+2)%3;
        int i2 = (i1+2)%3;
        a[0] -= F1(i0,0)*F1(i1,1)*F1(i2,2);
        a[1] -= F2(i0,0)*F1(i1,1)*F1(i2,2)+
                F1(i0,0)*F2(i1,1)*F1(i2,2)+
                F1(i0,0)*F1(i1,1)*F2(i2,2);
        a[2] -= F1(i0,0)*F2(i1,1)*F2(i2,2)+
                F2(i0,0)*F1(i1,1)*F2(i2,2)+
                F2(i0,0)*F2(i1,1)*F1(i2,2);
        a[3] -= F2(i0,0)*F2(i1,1)*F2(i2,2);
    }
}

// Compute the real roots of a third order polynomial.
// Return the number of roots found (1 or 3).
static int FindCubicRoots(float coeff[4], float x[3])
{
    float a1 = coeff[2] / coeff[3];
    float a2 = coeff[1] / coeff[3];
    float a3 = coeff[0] / coeff[3];

    double Q = (a1 * a1 - 3 * a2) / 9;
    double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
    double Q3 = Q * Q * Q;
    double d = Q3 - R * R;

    if (d >= 0) { // 3 real roots
        double theta = acos(R / sqrt(Q3));
        double sqrtQ = sqrt(Q);
        x[0] = -2 * sqrtQ * cos( theta             / 3) - a1 / 3;
        x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
        x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
        return 3;
    } else { // 1 real root
        double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
        if (R > 0)
            e = -e;
        x[0] = (e + Q / e) - a1 / 3.;
        return 1;
    }
}

// Epipolar geometry associated to 7 correspondence pairs.
// \param[in]  match The full set of correspondences
// \param[in]  k     Indices of the 7 correspondences
// \param[out] F     1 or 3 fundamental matrices
static void epipolar(const std::vector<Match>& match, const size_t k[7],
                     std::vector<libNumerics::matrix<float> >& F)
{
    libNumerics::matrix<double> c(7,9); // Build 7x9 matrix from point matches
    for(int i=0; i < 7; i++) {
        const Match& m = match[k[i]];
        c(i,0) = m.x1*m.x2;
        c(i,1) = m.y1*m.x2;
        c(i,2) =      m.x2;
        c(i,3) = m.x1*m.y2;
        c(i,4) = m.y1*m.y2;
        c(i,5) =      m.y2;
        c(i,6) = m.x1;
        c(i,7) = m.y1;
        c(i,8) = 1.0;
    }
    libNumerics::SVD svd(c);

    // Build basis of solutions
    libNumerics::matrix<float> F1(3,3), F2(3,3);
    F1.read( svd.V().col(7) );
    F2.read( svd.V().col(8) );
    F2 -= F1;

    float a[4]; // Build cubic polynomial P(x)=det(F1+xF2)
    characPoly(F1, F2, a);

    float z[3];
    F.clear();
    for(int m = FindCubicRoots(a,z)-1; m >=0; m--)
        F.push_back(F1+z[m]*F2);
}

// Store square distance to epipolar line and index of match
typedef std::pair<float,size_t> ErrorIndex;

// Increasing sequence of squared distances to epipolar lines with index
static void sortErrors(const libNumerics::matrix<float>& F,
                       const std::vector<Match>& m,
                       std::vector<ErrorIndex>& e)
{
    e.clear();
    size_t n = m.size();
    for(size_t i=0; i < n; i++) {
        float x = m[i].x1;
        float y = m[i].y1;
        float a = F(0,0)*x+F(0,1)*y+F(0,2);
        float b = F(1,0)*x+F(1,1)*y+F(1,2);
        float c = F(2,0)*x+F(2,1)*y+F(2,2);
        float d = (a*m[i].x2+b*m[i].y2+c);
        e.push_back( ErrorIndex((d*d)/(a*a+b*b),i) );
    }
    std::sort(e.begin(), e.end());
}

// Find best NFA and its index wrt square error threshold in e.
static ErrorIndex bestNFA(const std::vector<ErrorIndex>& e,
                          float logalpha0, float loge0,
                          const float* logcn, const float* logc7)
{
    const size_t n = e.size();
    ErrorIndex best(FLT_MAX, 6);
    for(size_t i=7; i < n; i++) {
        float logalpha = logalpha0+0.5f*(float)log10(e[i].first);
        ErrorIndex ei(loge0+logalpha*(float)(i-6)+logcn[i+1]+logc7[i+1], i);
        if(ei.first < best.first)
            best = ei;
    }
    return best;
}

/// ORSA algorithm.
libNumerics::matrix<float>
orsa(const std::vector<Match>& match,
     int t, bool verb, int mode, bool stop, float logalpha0,
     std::vector<size_t>& inliers, float& errorMax)
{
    // check size
    size_t n = match.size();
    if (n < 14) {
        std::cerr << "Inconsistent sizes." << std::endl;
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

    // tabulate logcombi
    float loge0 = log10(3.0f*(n-7));
    float* logcn = makelogcombi_n(n);
    float* logc7 = makelogcombi_k(7,n); 

    if(mode==3) // choose mode
        mode = ((logcn[7]<=log10((float)t))? 0: 2);
    if(verb)
        switch(mode) {
        case 0:
            std::cout << "Deterministic mode (systematic search). Test "
                      << int(0.5f+pow(10.f,logc7[n])) <<" bases" <<std::endl;
            break;
        case 1:
            std::cout << "Pure stochastic mode, no optimization." <<std::endl;
            break;
        case 2:
            std::cout << "Optimized stochastic mode (ORSA)." <<std::endl;
            break;
        default:
            assert(false);
        }

    std::vector<libNumerics::matrix<float> > F; // Up to 3 candidates per sample
    libNumerics::matrix<float> Fopti(3,3); // Optimal F matrix
    std::vector<ErrorIndex> e;
    e.reserve(n);

    std::vector<size_t> id(n);
    for(size_t i=0; i < n; i++)
        id[i] = i;

    float minNFA = FLT_MAX;
    int maxniter = (mode==2? t-t/10: t); // Reserve 10% for optimized sampling
    int niter = 0;

    size_t k[8]; // Sample of 7 indices+sentinel value
    if(mode==0) {
        for(size_t i=0; i < 7; i++)
            k[i]=i;
        k[6] = 5; // so that we will get k[6]=6 for first sample
        k[7] = n; // sentinel value
    }

    /********** MAIN LOOP **********/
    bool cont=true, optim=false;
    do {
      niter++;

      // Build new sample of 7 correspondences
      if(mode)
          random_p7(k, id.size());
      else if(! next_sample(k))
          break;

      // Find epipolar transform
      size_t idk[7];
      for(int i=0; i<7; i++)
          idk[i] = id[k[i]];
      epipolar(match, idk, F);

      std::vector<libNumerics::matrix<float> >::const_iterator it=F.begin();
      for(; it != F.end(); ++it) { // Loop on 1 or 3 solutions
          sortErrors(*it, match, e);
          ErrorIndex best = bestNFA(e, logalpha0, loge0, logcn, logc7);

          bool better = (best.first<minNFA);
          if(better) {
              minNFA = best.first;
              inliers.clear();
              for(size_t i=0; i <= best.second; i++)
                  inliers.push_back(e[i].second);
              errorMax = e[best.second].first; // Error threshold
              Fopti = *it;
          }

          // ORSA optimization: draw samples among best set of inliers so far
          if(mode==2 &&
             ((better && minNFA<0.0f) || (niter==maxniter && !optim))) {
              if(!optim) { // First time here: add back 10% iterations
                  optim = true;
                  maxniter = t;
              }
              id = inliers;
              if(verb)
                  std::cout << "   nfa=" << minNFA
                            << " size=" << id.size()
                            << " (niter=" << niter << ")" <<std::endl;
          }
      }

      if(stop && minNFA<0.0f)
          cont = false;
      else if(mode)
          cont = (niter<maxniter);
    } while(cont);

    if(verb)
        std::cout << "best matching found:  " << inliers.size()
                  << " points  log(nfa)=" << minNFA
                  << "  (" << niter << " iterations)" <<std::endl;

    delete [] logc7;
    delete [] logcn;
    errorMax = std::sqrt(errorMax);
    return Fopti;
}
