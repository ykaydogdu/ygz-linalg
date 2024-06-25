#include "ygzlinalg.hpp"
#include <stdexcept>

/* ************************************************************************* */
/*                        Constructors and Destructor                        */
/* ************************************************************************* */
// Default constructor
template <class T>
ygzVector<T>::ygzVector()
{
    nDims = 0;
    data = nullptr;
}

// Constructor with size
template <class T>
ygzVector<T>::ygzVector(size_t nDims)
{
    this->nDims = nDims;
    data = new T[nDims];
}

// Constructor with input data
template <class T>
ygzVector<T>::ygzVector(size_t nDims, const T *inputData)
{
    this->nDims = nDims;
    data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        data[i] = inputData[i];
    }
}

// Constructor with std::vector input data
template <class T>
ygzVector<T>::ygzVector(vector<T> &inputData)
{
    nDims = inputData.size();
    data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        data[i] = inputData.at(i);
    }
}

// Destructor
template <class T>
ygzVector<T>::~ygzVector()
{
    if (data != nullptr)
    {
        delete[] data;
    }
}

/* ************************************************************************* */
/*                             Property Access                               */
/* ************************************************************************* */
// Get the number of dimensions
template <class T>
size_t ygzVector<T>::getNumDims() const
{
    return nDims;
}

// Get the element at the specified index
template <class T>
T ygzVector<T>::getElement(int index) const
{
    if (index < 0 || index >= nDims)
    {
        throw std::invalid_argument("Index out of bounds");
    }
    return data[index];
}

// Set the element at the specified index to the given value
template <class T>
void ygzVector<T>::setElement(int index, T value)
{
    if (index < 0 || index >= nDims)
    {
        throw std::invalid_argument("Index out of bounds");
    }
    data[index] = value;
}

/* ************************************************************************* */
/*                              Computations                                 */
/* ************************************************************************* */
// Compute the norm of the vector
template <class T>
T ygzVector<T>::norm(int l) const
{
    if (l != 1 && l != 2)
        throw std::invalid_argument("Only L1 and L2 norms are supported");
    T result = 0;
    if (l == 1)
    {
        for (size_t i = 0; i < nDims; i++)
        {
            result += fabs(data[i]);
        }
        return result;
    } else 
    {
        for (size_t i = 0; i < nDims; i++)
        {
            result += data[i] * data[i];
        }
        return sqrt(result);
    }
}

// Return a normalized vector
template <class T>
ygzVector<T> ygzVector<T>::normalized() const
{
    T n = norm();
    if (n == 0)
    {
        throw std::invalid_argument("Cannot normalize a zero vector");
    }

    ygzVector<T> result(nDims, data);
    return result * (static_cast<T>(1.0) / n);
}

// Normalize the vector in place
template <class T>
void ygzVector<T>::normalize()
{
    T n = norm();
    if (n == 0)
    {
        throw std::invalid_argument("Cannot normalize a zero vector");
    }

    for (size_t i = 0; i < nDims; i++)
    {
        data[i] /= n;
    }
}

/* ************************************************************************* */
/*                                Operations                                 */
/* ************************************************************************* */
// Vector addition
template <class T>
ygzVector<T> ygzVector<T>::operator+(const ygzVector<T> &rhs) const
{
    if (nDims != rhs.nDims)
    {
        throw std::invalid_argument("Vector dimensions do not match");
    }

    ygzVector<T> result;
    result.nDims = nDims;
    result.data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        result.data[i] = data[i] + rhs.data[i];
    }

    return result;
}

// Vector subtraction
template <class T>
ygzVector<T> ygzVector<T>::operator-(const ygzVector<T> &rhs) const
{
    if (nDims != rhs.nDims)
    {
        throw std::invalid_argument("Vector dimensions do not match");
    }

    ygzVector<T> result;
    result.nDims = nDims;
    result.data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        result.data[i] = data[i] - rhs.data[i];
    }

    return result;
}

// Scalar multiplication
template <class T>
ygzVector<T> ygzVector<T>::operator*(const T &rhs) const
{
    ygzVector<T> result;
    result.nDims = nDims;
    result.data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        result.data[i] = data[i] * rhs;
    }

    return result;
}

// Scalar multiplication (friend function)
template <class T>
ygzVector<T> operator*(const T &lhs, const ygzVector<T> &rhs)
{
    return rhs * lhs;
}

// Dot product
template <class T>
T ygzVector<T>::dotProduct(const ygzVector<T> &lhs, const ygzVector<T> &rhs)
{
    if (lhs.nDims != rhs.nDims)
    {
        throw std::invalid_argument("Vector dimensions do not match");
    }

    T result = 0;
    for (size_t i = 0; i < lhs.nDims; i++)
    {
        result += lhs.data[i] * rhs.data[i];
    }

    return result;
}

// Cross product
template <class T>
ygzVector<T> ygzVector<T>::crossProduct(const ygzVector<T> &lhs, const ygzVector<T> &rhs)
{
    if (lhs.nDims != 3 || rhs.nDims != 3)
    {
        throw std::invalid_argument("Cross product is only defined for 3D vectors");
    }

    ygzVector<T> result;
    result.nDims = 3;
    result.data = new T[3];
    result.data[0] = lhs.data[1] * rhs.data[2] - lhs.data[2] * rhs.data[1];
    result.data[1] = lhs.data[2] * rhs.data[0] - lhs.data[0] * rhs.data[2];
    result.data[2] = lhs.data[0] * rhs.data[1] - lhs.data[1] * rhs.data[0];

    return result;
}

/* ************************************************************************* */
/*                             Function Mapping                              */
/* ************************************************************************* */
// Apply a function to each element of the vector
template <class T>
ygzVector<T> ygzVector<T>::map(T (*f)(T)) const
{
    ygzVector<T> result;
    result.nDims = nDims;
    result.data = new T[nDims];
    for (size_t i = 0; i < nDims; i++)
    {
        result.data[i] = f(data[i]);
    }

    return result;
}

// Apply a function to each element of the vector in place
template <class T>
void ygzVector<T>::mapInPlace(T (*f)(T))
{
    for (size_t i = 0; i < nDims; i++)
    {
        data[i] = f(data[i]);
    }
}


