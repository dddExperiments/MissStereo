#define _USE_MATH_DEFINES // For Windows
#include "libNumerics/numerics.h"
#include "libNumerics/homography.h"
#include "libMatch/match.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace libNumerics;

/// Iterative minimization of Sampson error by Levenberg-Marquardt
class LMRectify : public MinLM {
public:
    LMRectify(const std::vector<Match>& m, int w, int h)
    : m_matches(m), m_w(w), m_h(h) {}
    virtual void modelData(const vector<flnum>& P,
                           vector<flnum>& ymodel) const;
    virtual void modelJacobian(const vector<flnum>& P,
                               matrix<flnum>& J) const;

    flnum f(const vector<flnum>& P) const;
    matrix<flnum> getK(const vector<flnum>& P) const;
    matrix<flnum> getInvK(const vector<flnum>& P) const;
    matrix<flnum> getInvKd(const vector<flnum>& P) const;
    matrix<flnum> getF(const vector<flnum>& P) const;
    matrix<flnum> getRl(const vector<flnum>& P) const;
    matrix<flnum> getRr(const vector<flnum>& P) const;

    void fix_center(Homography& Hl, Homography& Hr, const vector<flnum>&) const;
private:
    const std::vector<Match>& m_matches; ///< Set of matches
    int m_w, m_h; ///< Width/height of image
};

/// Focal length, encoded in P(5)
inline flnum LMRectify::f(const vector<flnum>& P) const
{ return std::pow((flnum)3.0, P(5)) * (m_w+m_h); }

/// K
matrix<flnum> LMRectify::getK(const vector<flnum>& P) const
{
    flnum f = this->f(P);
    matrix<double> K(3,3);
    K = 0.0;
    K(0,0) = K(1,1) = f;
    K(0,2) = .5*m_w;
    K(1,2) = .5*m_h;
    K(2,2) = 1.0;
    return K;
}

/// K^{-1}
matrix<flnum> LMRectify::getInvK(const vector<flnum>& P) const
{
    flnum f = this->f(P);
    matrix<flnum> K(3,3);
    K = 0;
    K(0,0) = K(1,1) = 1/f;
    K(2,2) = 1;
    K(0,2) = -m_w / (2*f);
    K(1,2) = -m_h / (2*f);
    return K;
}

/// Derivative of K^{-1}
matrix<flnum> LMRectify::getInvKd(const vector<flnum>& P) const
{
    flnum f = this->f(P);
    static flnum log3 = std::log((flnum)3.0);
    matrix<flnum> K(3,3);
    K = 0;
    K(0,0) = K(1,1) = -log3/f;
    K(0,2) = +m_w*log3 / (2*f);
    K(1,2) = +m_h*log3 / (2*f);
    return K;
}

/// Get rotation matrix of angle \a theta around \a axis ('x','y' or 'z').
/// \a deriv indicates whether we want the derivative of rotation matrix.
static matrix<flnum> getR(flnum theta, char axis, bool deriv=false)
{
    matrix<flnum> R(3,3);
    R=0;
    if(deriv)
        theta += .5*M_PI;
    flnum c=std::cos(theta), s=std::sin(theta);
    switch(axis) {
    case 'x':
        if(! deriv) R(0,0) = (flnum)1;
        R(1,1) = R(2,2) = c;
        R(1,2) = -(R(2,1) = s);
        break;
    case 'y':
        if(! deriv) R(1,1) = (flnum)1;
        R(0,0) = R(2,2) = c;
        R(0,2) = -(R(2,0) = s);
        break;
    case 'z':
        if(! deriv) R(2,2) = (flnum)1;
        R(0,0) = R(1,1) = c;
        R(0,1) = -(R(1,0) = s);
        break;
    default: assert(false);
    }
    return R;
}

/// Left rotation matrix of rectification
matrix<flnum> LMRectify::getRl(const vector<flnum>& P) const
{
    flnum ay=P(0), az=P(1);
    matrix<flnum> Ry=getR(ay,'y');
    matrix<flnum> Rz=getR(az,'z');
    return Rz*Ry;
}

/// Right rotation matrix of rectification
matrix<flnum> LMRectify::getRr(const vector<flnum>& P) const
{
    flnum ax=P(2), ay=P(3), az=P(4);
    matrix<flnum> Rx=getR(ax,'x');
    matrix<flnum> Ry=getR(ay,'y');
    matrix<flnum> Rz=getR(az,'z');
    return Rz*Ry*Rx;
}

/// Get fundamental matrix
matrix<flnum> LMRectify::getF(const vector<flnum>& P) const
{
    matrix<flnum> invK=getInvK(P);
    matrix<flnum> F(3,3);
    F = 0; F(1,2) = -(F(2,1) = 1);
    return (getRl(P)*invK).t() * F * (getRr(P)*invK);
}

/// Sampson error
void LMRectify::modelData(const vector<flnum>& P,
                          vector<flnum>& ymodel) const
{
    matrix<flnum> invK=getInvK(P);
    matrix<flnum> F(3,3);
    F = 0; F(1,2) = -(F(2,1) = 1);
    F = (getRl(P)*invK).t() * F * (getRr(P)*invK);

    std::vector<Match>::const_iterator it = m_matches.begin();
    for(int i=0; it != m_matches.end(); ++it, ++i) {
        vector<flnum> ml(3), mr(3);
        ml(0) = it->x1; ml(1) = it->y1; ml(2) = 1;
        mr(0) = it->x2; mr(1) = it->y2; mr(2) = 1;
        ymodel(i) =
            dot(ml,F*mr) /
            std::sqrt((F*mr).copy(0,1).qnorm() + (F.t()*ml).copy(0,1).qnorm());
    }
}

/// Jacobian of sampson error
void LMRectify::modelJacobian(const libNumerics::vector<flnum>& P,
                              libNumerics::matrix<flnum>& J) const
{
    const matrix<flnum> invK=getInvK(P), dinvK=getInvKd(P);
    const matrix<flnum> Rly=getR(P(0),'y');
    const matrix<flnum> Rlz=getR(P(1),'z');
    const matrix<flnum> Rrx=getR(P(2),'x');
    const matrix<flnum> Rry=getR(P(3),'y');
    const matrix<flnum> Rrz=getR(P(4),'z');
    matrix<flnum> F(3,3), Fl(3,3), Fr(3,3);
    F = 0; F(1,2) = -(F(2,1) = 1);
    Fl = (getRl(P)*invK).t() * F;
    Fr = F * (getRr(P)*invK);
    F = (getRl(P)*invK).t() * F * (getRr(P)*invK);

    const matrix<flnum> dF0 = (Rlz*getR(P(0),'y',true)*invK).t()*Fr;
    const matrix<flnum> dF1 = (getR(P(1),'z',true)*Rly*invK).t()*Fr;
    const matrix<flnum> dF2 = Fl*Rrz*Rry*getR(P(2),'x',true)*invK;
    const matrix<flnum> dF3 = Fl*Rrz*getR(P(3),'y',true)*Rrx*invK;
    const matrix<flnum> dF4 = Fl*getR(P(4),'z',true)*Rry*Rrx*invK;
    const matrix<flnum> dF5 = (getRl(P)*dinvK).t() * Fr + Fl * (getRr(P)*dinvK);

    std::vector<Match>::const_iterator it = m_matches.begin();
    for(int i=0; it != m_matches.end(); ++it, ++i) {
        vector<flnum> ml(3), mr(3);
        ml(0) = it->x1; ml(1) = it->y1; ml(2) = 1;
        mr(0) = it->x2; mr(1) = it->y2; mr(2) = 1;
        flnum num=dot(ml,F*mr);
        vector<flnum> fmr=(F*mr).copy(0,1), ftml=(F.t()*ml).copy(0,1);
        flnum denom=std::sqrt(fmr.qnorm()+ftml.qnorm());

#define DERIV(dF)                   \
        dot(ml,dF*mr)/denom - \
        num*(dot(fmr,(dF*mr).copy(0,1))+ \
             dot(ftml,(dF.t()*ml).copy(0,1)))   \
        /(denom*denom*denom)

        J(i,0) = DERIV(dF0);
        J(i,1) = DERIV(dF1);
        J(i,2) = DERIV(dF2);
        J(i,3) = DERIV(dF3);
        J(i,4) = DERIV(dF4);
        J(i,5) = DERIV(dF5);
#undef DERIV
    }
}

/// Change x-angle and principal points in \a Hl and \a Hr so as to fix the
/// point (w/2,h/2).
void LMRectify::fix_center(Homography& Hl, Homography& Hr,
                           const vector<flnum>& P) const
{
    matrix<double> K=getK(P), invK=getInvK(P);

    // Rotate around epipolar axis to keep center point at same ordinate
    double x=.5*m_w, y=.5*m_h;
    Hl(x,y);
    double alpha=std::atan2(y-.5*m_h,this->f(P));
    Hl.mat() = K*getR(alpha,'x')*getRl(P)*invK;
    Hr.mat() = K*getR(alpha,'x')*getRr(P)*invK;

    // Move principal point to keep center point at same abscissa in each image
    x=.5*m_w, y=.5*m_h;
    Hl(x,y);
    K(0,2)=m_w-x;
    Hl.mat() = K*getR(alpha,'x')*getRl(P)*invK;
    x=.5*m_w, y=.5*m_h;
    Hr(x,y);
    K(0,2)=m_w-x;
    Hr.mat() = K*getR(alpha,'x')*getRr(P)*invK;
}

/// Compute rectifying homographies.
/// Rectification invariant degrees of freedom are computed so as to keep image
/// centers fixed.
std::pair<float,float> compRectif(int w, int h, const std::vector<Match>& m,
                                  Homography& Hl, Homography& Hr)
{
    LMRectify lm(m, w, h);
    vector<flnum> P(6), ydata((int)m.size());
    P = 0;
    std::pair<float,float> err;
    lm.modelData(P, ydata);
    err.first = static_cast<float>( sqrt(ydata.qnorm()/ydata.nrow()) );
    ydata = 0;
    err.second = (float)lm.minimize(P, ydata);
    matrix<double> K = lm.getK(P), invK = lm.getInvK(P);
    Hl.mat() = K*lm.getRl(P)*invK;
    Hr.mat() = K*lm.getRr(P)*invK;

    lm.fix_center(Hl, Hr, P);

    std::cout <<"LM iterations: " <<lm.iterations <<" f=" <<K(0,0) <<std::endl;
    return err;
}

/// Compute and print min and max disparity
void printDisparity(const std::vector<Match>& match,
                    const libNumerics::Homography& Hl,
                    const libNumerics::Homography& Hr)
{
    std::vector<Match>::const_iterator it=match.begin();
    double min=DBL_MAX, max=-DBL_MAX;
    for(; it != match.end(); ++it) {
        double xl=it->x1, yl=it->y1;
        Hl(xl,yl);
        double xr=it->x2, yr=it->y2;
        Hr(xr,yr);
        xr -= xl;
        if(xr < min)
            min = xr;
        if(xr > max)
            max = xr;
    }
    std::cout << "Disparity: "
              << (int)floor(min) << " " << (int)ceil(max) << std::endl;
}

int main(int argc, char** argv)
{
    if(argc != 6) {
        std::cerr << "Usage: " << argv[0] << " match.txt w h Hl Hr" <<std::endl;
        return 1;
    }

    std::vector<Match> match;
    if(! loadMatch(argv[1],match)) {
        std::cerr << "Failed reading " << argv[1] << std::endl;
        return 1;
    }

    int w=0,h=0;
    if(! (std::istringstream(argv[2]) >> w).eof()) w=0;
    if(! (std::istringstream(argv[3]) >> h).eof()) h=0;
    if(w <=0 || h <= 0) {
        std::cerr << "Wrong dimensions of image" << std::endl;
        return 1;
    }

    libNumerics::Homography Hl, Hr;
    std::pair<float,float> e = compRectif(w, h, match, Hl, Hr);
    std::cout << "Initial rectification error: " <<e.first <<" pix" <<std::endl;
    std::cout << "Final rectification error: " << e.second <<" pix" <<std::endl;
    printDisparity(match, Hl, Hr);

    // Output files
    std::ofstream f1(argv[4]), f2(argv[5]);
    if((f1 << Hl.mat() << std::endl).fail()) {
        std::cerr << "Error writing file " << argv[4] << std::endl;
        return 1;
    }
    if((f2 << Hr.mat() << std::endl).fail()) {
        std::cerr << "Error writing file " << argv[5] << std::endl;
        return 1;
    }

    return 0;
}
