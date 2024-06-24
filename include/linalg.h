#ifndef LINALG_H
#define LINALG_H

template <class T>
class ygzMatrix
{
public:
    // Constructors
    ygzMatrix();
    ygzMatrix(size_t nRows, size_t nCols);
    ygzMatrix(size_t nRows, size_t nCols, const T* val);
    ygzMatrix(const ygzMatrix<T> &copyMatrix); // Copy constructor
    
    // Destructor
    ~ygzMatrix();

    // Resizing
    bool resize(size_t nRows, size_t nCols);

    // Element access
    T getElement(size_t nRow, size_t nCol) const;
    bool setElement(size_t nRow, size_t nCol, T value);
    int getNumRows() const;
    int getNumCols() const;

    // Equality
    bool operator==(const ygzMatrix<T> &rhs) const;

    // Arithmetic operations
    template <class U> friend ygzMatrix<U> operator+(const ygzMatrix<U> &lhs, const ygzMatrix<U> &rhs);
    template <class U> friend ygzMatrix<U> operator+(const ygzMatrix<U> &matrix, const U &scalar);
    template <class U> friend ygzMatrix<U> operator+(const U &scalar, const ygzMatrix<U> &rhs);

    template <class U> friend ygzMatrix<U> operator-(const ygzMatrix<U> &lhs, const ygzMatrix<U> &rhs);
    template <class U> friend ygzMatrix<U> operator-(const ygzMatrix<U> &matrix, const U &scalar);
    template <class U> friend ygzMatrix<U> operator-(const U &scalar, const ygzMatrix<U> &matrix);

    template <class U> friend ygzMatrix<U> operator*(const ygzMatrix<U> &lhs, const ygzMatrix<U> &rhs);
    template <class U> friend ygzMatrix<U> operator*(const ygzMatrix<U> &matrix, const U &scalar);
    template <class U> friend ygzMatrix<U> operator*(const U &scalar, const ygzMatrix<U> &matrix);

    // Special matrix generation
    static ygzMatrix random(size_t nRows, size_t nCols);
    static ygzMatrix identity(size_t n);

private:
    int sub2Ind(int r, int c) const;

private:
    int nRows, nCols, nElements;
    T *data;
};

#endif