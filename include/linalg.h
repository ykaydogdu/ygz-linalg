#ifndef LINALG_H
#define LINALG_H

#include <cstddef>
#include <vector>
using std::vector;

template <class T>
class ygzMatrix
{
public:
    // Constructors
    ygzMatrix();
    ygzMatrix(size_t nRows, size_t nCols);
    ygzMatrix(size_t nRows, size_t nCols, const T* inputData);
    ygzMatrix(size_t nRows, size_t nCols, const vector<T> *inputData);
    ygzMatrix(const ygzMatrix<T> &copyMatrix); // Copy constructor
    
    // Destructor
    ~ygzMatrix();

    // Configuration
    bool resize(size_t nRows, size_t nCols);
    void setToIdentity();

    // Element access
    T getElement(size_t nRow, size_t nCol) const;
    bool setElement(size_t nRow, size_t nCol, T value);
    int getNumRows() const;
    int getNumCols() const;

    // Operations
    bool inverseInPlace();

    // Equality
    bool operator==(const ygzMatrix<T> &rhs) const;
    bool compare(const ygzMatrix<T> & m2, double tolerance) const;

    // Arithmetic operations
    template <class U> friend ygzMatrix<U> operator+(const ygzMatrix<U> &lhs, const ygzMatrix<U> &rhs);
    template <class U> friend ygzMatrix<U> operator+(const ygzMatrix<U> &matrix, const U &scalar);
    template <class U> friend ygzMatrix<U> operator+(const U &scalar, const ygzMatrix<U> &matrix);

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
    bool isSquare() const;
    bool closeEnough(T a, T b, double tolerance) const;
    void swapRow(int i, int j);
    void multRow(int i, T factor);
    void multAddRow(int i, int j, T factor);
    // void swapCol(int i, int j);
    // void multCol(int i, T factor);
    // void multAddCol(int i, int j, T factor);
    ygzMatrix<T> join(const ygzMatrix<T> &m);
    void separate(ygzMatrix<T> *m1, ygzMatrix<T> *m2, int colNum) const;
    int findRowWithMax(int col, int startingRow) const; // pivot finding

private:
    int nRows, nCols, nElements;
    T *data;
};

#endif