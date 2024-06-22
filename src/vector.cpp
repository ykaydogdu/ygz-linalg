#include "vector.h"
#include "matrix.h"
#include <math.h>
#include <stdexcept>

// default constructor
Vector::Vector()
{
    n = 0;
    val = nullptr;
}

// initialize a column vector with given dimensions
Vector::Vector(int n)
{
    this->n = n;
    val = new double[n]; 
}

// initialize a column vector from a given 1D array
Vector::Vector(double *val)
{
    this->n = sizeof(val) / sizeof(val[0]);
    this->val = new double[n];
    for (int i = 0; i < n; i++)
    {
        this->val[i] = val[i];
    }
}

// destructor
Vector::~Vector()
{
    delete[] val;
}

// implicit conversion from matrix to vector
Matrix::operator Vector()
{
    if (c != 1)
    {
        throw std::invalid_argument("Cannot convert a non-column matrix to a vector");
    }
    double* res = new double[r];
    for (int i = 0; i < r; i++)
    {
        res[i] = val[i][0];
    }
    return Vector(res);
}

Vector operator*(double const scalar, Vector const &m)
{
    Vector result = Vector(m.n);
    for (int i = 0; i < m.n; i++)
    {
        result.val[i] = m.val[i] * scalar;
    }
    return result;
}

Vector::operator Matrix()
{
    double** res = new double*[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = new double[1];
        res[i][0] = val[i];
    }
    return Matrix(res);
}

double Vector::norm(int l)
{
    double result = 0;
    if (l == 1)
    {
        for (int i = 0; i < n; i++)
        {
            result += abs(val[i]);
        }
        return result;
    }
    else if (l == 2)
    {
        for (int i = 0; i < n; i++)
        {
            result += val[i] * val[i];
        }
        return sqrt(result);
    }
    else
    {
        throw std::invalid_argument("Invalid norm type");
    }
}

Vector Vector::normalize()
{
    double n = norm();
    if (n == 0)
    {
        throw std::invalid_argument("Cannot normalize a zero vector");
    }
    return (1 / n) * *this;
}