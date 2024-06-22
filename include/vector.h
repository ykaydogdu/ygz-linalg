#ifndef VECTOR_H
#define VECTOR_H

#include "matrix.h"
#include <math.h>

class Vector {
public:
    Vector(int n);
    Vector(double* val);
    ~Vector();
    friend Vector operator*(Vector const& m1, Vector const& m2);
    friend Vector operator*(double const scalar, Vector const& m);
    operator Matrix();
    double norm(int l = 2);
    Vector normalize();
private:
    int n;
    double* val;
};

#endif