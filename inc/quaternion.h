#ifndef QUATERNION_H
#define QUATERNION_H
#include<array>

struct quat {
    
    double r;
    double i;
    double j;
    double k;
    
    inline quat& operator+ (const quat&);
    inline quat& operator* (const quat&);
    inline quat& operator/ (const double);
    inline quat& operator* (const double);
    inline quat inv (void) const;
    inline quat conj (void) const;
    inline double mag (void) const;
    
};

int rotate(std::array<double,3>& r,const std::array<double,3>& axis, const double angle);


#endif