#include "matrix.h"
#include <stdexcept>

// initialize a matrix with given dimensions
Matrix::Matrix(int r, int c)
{
    this->r = r;
    this->c = c;
    val = new double*[r];
    for (int i = 0; i < r; i++)
    {
        val[i] = new double[c];
    }
}

// initialize a matrix from a given 2D array
Matrix::Matrix(double** val)
{
    this->r = sizeof(val) / sizeof(val[0]);
    this->c = sizeof(val[0]) / sizeof(val[0][0]);
    this->val = new double*[r];
    for (int i = 0; i < r; i++)
    {
        this->val[i] = new double[c];
        for (int j = 0; j < c; j++)
        {
            this->val[i][j] = val[i][j];
        }
    }
}

// generate a random matrix with entries between [0, 1] in given dimensions
Matrix Matrix::random(int r, int c)
{
    Matrix result(r, c);
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            result.val[i][j] = rand() / (double)RAND_MAX;
        }
    }
    return result;
}

// generate an identity matrix of given size
Matrix Matrix::identity(int n)
{
    Matrix result(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            result.val[i][j] = (i == j) ? 1 : 0;
        }
    }
    return result;
}

Matrix::~Matrix()
{
    delete(val);
}

Matrix operator+(Matrix const& m1, Matrix const& m2)
{
    // check if the matrices have the same dimensions
    if (m1.r != m2.r || m1.c != m2.c)
    {
        throw std::invalid_argument("The matrices must have the same dimensions");
    }

    Matrix result(m1.r, m1.c);
    for (int i = 0; i < m1.r; i++)
    {
        for (int j = 0; j < m1.c; j++)
        {
            result.val[i][j] = m1.val[i][j] + m2.val[i][j];
        }
    }
    return result;
}

// performs the matrix multiplication
Matrix operator*(Matrix const& m1, Matrix const& m2)
{   
    // check if the matrices can be multiplied
    if (m1.c != m2.r)
    {
        throw std::invalid_argument("The matrices cannot be multiplied");
    }

    Matrix result(m1.r, m2.c);
    for (int i = 0; i < m1.r; i++)
    {
        for (int j = 0; j < m2.c; j++)
        {
            result.val[i][j] = 0;
            for (int k = 0; k < m1.c; k++)
            {
                result.val[i][j] += m1.val[i][k] * m2.val[k][j];
            }
        }
    }
    return result;
}

// performs scalar multiplication
Matrix operator*(double scalar, Matrix const& m)
{
    Matrix result(m.r, m.c);
    for (int i = 0; i < m.r; i++)
    {
        for (int j = 0; j < m.c; j++)
        {
            result.val[i][j] = m.val[i][j] * scalar;
        }
    }
    return result;
}

// returns the transpose of the matrix
Matrix Matrix::transpose()
{
    Matrix result(c, r);
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            result.val[j][i] = val[i][j];
        }
    }
    return result;
}