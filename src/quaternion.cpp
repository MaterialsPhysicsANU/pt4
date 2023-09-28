#include "quaternion.h"
#include "vec.h"

#include <cmath>

inline quat quat::inv(void) const{
    quat inv = *this;
    inv = inv.conj()/inv.mag();
    return inv;
}

inline quat quat::conj(void) const{
    return {r,-i,-j,-k};
}

inline double quat::mag(void) const{
    quat q = *this;
    double mag = sqrt(q.r*q.r + q.i*q.i + q.j*q.j + q.k*q.k);
    return mag;
}

inline quat& quat::operator* (const double a){
    r *= a;    
    i *= a;    
    j *= a;    
    k *= a;
    return *this;
}

inline quat& quat::operator/ (const double a){
    r /= a;    
    i /= a;    
    j /= a;    
    k /= a;
    return *this;
}

inline quat& quat::operator+ (const quat& q){
    r += q.r;    
    i += q.i;    
    j += q.j;    
    k += q.k;
    return *this;
}
    
inline quat& quat::operator* (const quat& q){
    double rt = r*q.r - i*q.i - j*q.j - k*q.k;    
    double it = r*q.i + i*q.r - k*q.j + j*q.k;    
    double jt = r*q.j + j*q.r + k*q.i - i*q.k;    
    double kt = r*q.k + k*q.r - j*q.i + i*q.j;
    r = rt;
    i = it;
    j = jt;
    k = kt;
    return *this;
}

int rotate(std::array<double,3>& r,const std::array<double,3>& axis, const double angle){
    double mg = mag(axis);
    if (mg == 0){return 1;}
    double mgi = 1.0/mg;
    std::array<double,3> uv = axis*mgi;
    double half_sin = sin(angle*0.5);
    double half_cos = cos(angle*0.5);
    quat qr = {0.0, r[0],r[1],r[2]}; 
    quat qu = {half_cos, half_sin*uv[0],half_sin*uv[1],half_sin*uv[2]}; 
    quat qui = qu.conj();//since qu is unitary
    qr = qu*qr*qui;
    
    //store results
    r[0] = qr.i;
    r[1] = qr.j;
    r[2] = qr.k;
    return 0;
}