#ifndef MATRIX_H
#define MATRIX_H

class Vector;

class Matrix {
public:
    Matrix(int r, int c);
    Matrix(double** val);
    static Matrix random(int r, int c);
    static Matrix identity(int n);
    ~Matrix();
    friend Matrix operator+(Matrix const& m1, Matrix const& m2);
    friend Matrix operator*(Matrix const& m1, Matrix const& m2);
    friend Matrix operator*(double scalar, Matrix const& m);
    Matrix transpose();
    operator Vector();
protected:
    int r, c;
    double** val;
};

#endif