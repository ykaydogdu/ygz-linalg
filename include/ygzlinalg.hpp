#ifndef LINALG_H
#define LINALG_H

#include <cstddef>
#include <vector>
using std::vector;

#define TOLERANCE 1e-5

template <class T>
class ygzVector
{
public:
    // Constructors
    ygzVector();
    ygzVector(size_t nDims);
    ygzVector(size_t nDims, const T* inputData);
    ygzVector(vector<T> &inputData);
    ~ygzVector();

    // getter
    size_t getNumDims() const;

    // element access
    T getElement(int index) const;
    void setElement(int index, T value);

    // Computations
    T norm(int l = 2) const;
    ygzVector<T> normalized() const; // return a normalized vector
    void normalize(); // normalize the vector in place

    // Operators
    ygzVector<T> operator+(const ygzVector<T> &rhs) const;
    ygzVector<T> operator-(const ygzVector<T> &rhs) const;
    ygzVector<T> operator*(const T &rhs) const;
    template <class U> friend ygzVector<U> operator*(const U &lhs, const ygzVector<U> &rhs);

    // Function mapping
    ygzVector<T> map(T (*f)(T)) const;
    void mapInPlace(T (*f)(T));

    // Static Functions
    static T dotProduct(const ygzVector<T> &lhs, const ygzVector<T> &rhs);
    static ygzVector<T> crossProduct(const ygzVector<T> &lhs, const ygzVector<T> &rhs);

private:
    int nDims;
    T *data;
};

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
    ygzMatrix<T> inverse() const;
    T determinant() const;

    // Equality
    bool operator==(const ygzMatrix<T> &rhs) const;
    bool operator!=(const ygzMatrix<T> &rhs) const;
    bool compare(const ygzMatrix<T> & m2, double tolerance = TOLERANCE) const;

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

    // matrix * vector
    template <class U> friend ygzVector<U> operator*(const ygzMatrix<U> &matrix, const ygzVector<U> &vector);
    template <class U> friend ygzMatrix<U> operator*(const ygzVector<U> &vector, const ygzMatrix<U> &r_vector);

    // Special matrix generation
    static ygzMatrix random(size_t nRows, size_t nCols);
    static ygzMatrix identity(size_t n);

private:
    int sub2Ind(int r, int c) const;
    bool isSquare() const;
    bool closeEnough(T a, T b, double tolerance = TOLERANCE) const;

    // Row and column operations
    void swapRow(int i, int j);
    void multRow(int i, T factor);
    void multAddRow(int i, int j, T factor);
    // void swapCol(int i, int j);
    // void multCol(int i, T factor);
    // void multAddCol(int i, int j, T factor);

    // Augmented matrix operations
    ygzMatrix<T> join(const ygzMatrix<T> &m);
    void separate(ygzMatrix<T> *m1, ygzMatrix<T> *m2, int colNum) const;

    int findRowWithMax(int col, int startingRow) const; // pivot finding

    ygzMatrix<T> findMinor(int row, int col) const; // find minor matrix

private:
    int nRows, nCols, nElements;
    T *data;
};

#endif