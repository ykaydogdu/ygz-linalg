#include "ygzlinalg.hpp"
#include <stdexcept>
#include <time.h>
#include <cmath>

#include <iostream>
template <class T>
void printMatrix(ygzMatrix<T> &m)
{
    for (int i = 0; i < m.getNumRows(); i++)
    {
        for (int j = 0; j < m.getNumCols(); j++)
        {
            std::cout << m.getElement(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/* ********************************************************************************************* */
/* 					    		Constructors & Destructor 										 */
/* ********************************************************************************************* */
// Default constructor
template <class T>
ygzMatrix<T>::ygzMatrix()
{
    nRows = 0;
    nCols = 0;
    nElements = 0;
    data = nullptr;
}

// Construct empty matrix (all elements are 0)
template <class T>
ygzMatrix<T>::ygzMatrix(size_t nRows, size_t nCols)
{
    this->nRows = nRows;
    this->nCols = nCols;
    nElements = nRows * nCols;
    data = new T[nElements];
    for (int i = 0; i < nElements; i++)
        data[i] = 0.0; // Initialize elements with default value
}

// Construct matrix from array
template <class T>
ygzMatrix<T>::ygzMatrix(size_t nRows, size_t nCols, const T* inputData)
{
    this->nRows = nRows;
    this->nCols = nCols;
    nElements = nRows * nCols;
    data = new T[nElements];
    for (int i = 0; i < nElements; i++)
        data[i] = inputData[i];
}

// Construct matrix from vector
template <class T>
ygzMatrix<T>::ygzMatrix(size_t nRows, size_t nCols, const vector<T> *inputData)
{
    this->nRows = nRows;
    this->nCols = nCols;
    nElements = nRows * nCols;
    data = new T[nElements];
    for (int i = 0; i < nElements; i++)
        data[i] = inputData->at(i);
}

// Copy constructor
template <class T>
ygzMatrix<T>::ygzMatrix(const ygzMatrix<T> &copyMatrix)
{
    nRows = copyMatrix.nRows;
    nCols = copyMatrix.nCols;
    nElements = nRows * nCols;
    data = new T[nElements];
    for (int i = 0; i < nElements; i++)
        data[i] = copyMatrix.data[i];
}

// Destructor
template <class T>
ygzMatrix<T>::~ygzMatrix()
{
    delete[] data;
}

/* ********************************************************************************************* */
/* 					    			  Special Matrix Generation									 */
/* ********************************************************************************************* */
// Generate a random matrix with values between 0 and 1
template <class T>
ygzMatrix<T> ygzMatrix<T>::random(size_t nRows, size_t nCols)
{
    srand(time(0)); // Seed the random number generator
    T* data = new T[nRows * nCols];
    for (int i = 0; i < nRows * nCols; i++)
        data[i] = static_cast<T>(rand()) / RAND_MAX; // Random value between 0 and 1

    ygzMatrix<T> randomMatrix(nRows, nCols, data);
    delete[] data;
    return randomMatrix;
}

// Generate an identity matrix
template <class T>
ygzMatrix<T> ygzMatrix<T>::identity(size_t n)
{
    T* data = new T[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            data[i * n + j] = i == j ? 1.0 : 0.0;
    }

    ygzMatrix<T> identityMatrix(n, n, data);
    delete[] data;
    return identityMatrix;
}

/* ********************************************************************************************* */
/* 					    		        Configuration 											 */
/* ********************************************************************************************* */
// Resize matrix
template <class T>
bool ygzMatrix<T>::resize(size_t nRows, size_t nCols)
{
    if (nRows == this->nRows && nCols == this->nCols)
        return true;

    delete[] data;
    this->nRows = nRows;
    this->nCols = nCols;
    nElements = nRows * nCols;
    data = new T[nElements];
    if (data == nullptr)
        return false;
    for (int i = 0; i < nElements; i++)
        data[i] = 0.0; // Initialize elements with default value

    return true;
}

// Set matrix to identity matrix
template <class T>
void ygzMatrix<T>::setToIdentity()
{
    if (!isSquare())
        throw std::invalid_argument("Matrix must be square to set to identity");

    for (int row = 0; row < nRows; row++)
    {
        for (int col = 0; col < nCols; col++)
        {
            setElement(row, col, row == col ? 1.0 : 0.0);
        }
    }
}

/* ********************************************************************************************* */
/* 					    			  Subscript Converter      									 */
/* ********************************************************************************************* */
// Convert subscripts to linear index
template <class T>
int ygzMatrix<T>::sub2Ind(int r, int c) const
{
    if (r < 0 || r >= nRows || c < 0 || c >= nCols)
        return -1;
    return r * nCols + c;
}

/* ********************************************************************************************* */
/* 					    			  Element access 											 */
/* ********************************************************************************************* */
// Get element at position (nRow, nCol)
template <class T>
T ygzMatrix<T>::getElement(size_t nRow, size_t nCol) const
{
    if (nRow >= nRows || nCol >= nCols)
        throw std::out_of_range("Index out of range");

    return data[sub2Ind(nRow, nCol)];
}

// Set element at position (nRow, nCol)
template <class T>
bool ygzMatrix<T>::setElement(size_t nRow, size_t nCol, T value)
{
    if (nRow >= nRows || nCol >= nCols)
        return false;

    data[sub2Ind(nRow, nCol)] = value;
    return true;
}

// Get number of rows
template <class T>
int ygzMatrix<T>::getNumRows() const
{
    return nRows;
}

// Get number of columns
template <class T>
int ygzMatrix<T>::getNumCols() const
{
    return nCols;
}

/* ********************************************************************************************* */
/* 					    			  Equality Operations										 */
/* ********************************************************************************************* */
// Check if two matrices are equal
template <class T>
bool ygzMatrix<T>::operator==(const ygzMatrix<T> &rhs) const
{
    if (nRows != rhs.nRows || nCols != rhs.nCols)
        return false;

    for (int i = 0; i < nElements; i++)
    {
        // if (data[i] != rhs.data[i])
        if (!closeEnough(data[i], rhs.data[i]))
            return false;
    }

    return true;
}

// Check if two matrices are not equal
template <class T>
bool ygzMatrix<T>::operator!=(const ygzMatrix<T> &rhs) const
{
    return !(*this == rhs);
}

// Compare two matrices with tolerance using MSE
template <class T>
bool ygzMatrix<T>::compare(const ygzMatrix<T> & m2, double tolerance) const
{
    if (nRows != m2.nRows || nCols != m2.nCols)
        return false;

    double sum = 0.0;
    for (int i = 0; i < nElements; i++)
    {
        sum += (data[i] - m2.data[i]) * (data[i] - m2.data[i]);
    }

    return sqrt(sum / (nElements - 1)) < tolerance;
}

/* ********************************************************************************************* */
/* 					    			  Arithmetic Operations										 */
/* ********************************************************************************************* */
// Matrix addition
template <class T>
ygzMatrix<T> operator+(const ygzMatrix<T> &lhs, const ygzMatrix<T> &rhs)
{
    if (lhs.nRows != rhs.nRows || lhs.nCols != rhs.nCols)
        throw std::invalid_argument("Matrices must have the same dimensions");

    T* result = new T[lhs.nElements];
    for (int i = 0; i < lhs.nElements; i++)
        result[i] = lhs.data[i] + rhs.data[i];

    ygzMatrix<T> resultMatrix(lhs.nRows, lhs.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Matrix-scalar addition
template <class T>
ygzMatrix<T> operator+(const ygzMatrix<T> &matrix, const T &scalar)
{
    T* result = new T[matrix.nElements];
    for (int i = 0; i < matrix.nElements; i++)
        result[i] = matrix.data[i] + scalar;

    ygzMatrix<T> resultMatrix(matrix.nRows, matrix.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Scalar-matrix addition
template <class T>
ygzMatrix<T> operator+(const T &scalar, const ygzMatrix<T> &matrix)
{
    return matrix + scalar; // Addition is commutative
}

// Matrix subtraction
template <class T>
ygzMatrix<T> operator-(const ygzMatrix<T> &lhs, const ygzMatrix<T> &rhs)
{
    if (lhs.nRows != rhs.nRows || lhs.nCols != rhs.nCols)
        throw std::invalid_argument("Matrices must have the same dimensions");

    T* result = new T[lhs.nElements];
    for (int i = 0; i < lhs.nElements; i++)
        result[i] = lhs.data[i] - rhs.data[i];

    ygzMatrix<T> resultMatrix(lhs.nRows, lhs.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Matrix-scalar subtraction
template <class T>
ygzMatrix<T> operator-(const ygzMatrix<T> &matrix, const T &scalar)
{
    T* result = new T[matrix.nElements];
    for (int i = 0; i < matrix.nElements; i++)
        result[i] = matrix.data[i] - scalar;

    ygzMatrix<T> resultMatrix(matrix.nRows, matrix.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Scalar-matrix subtraction
template <class T>
ygzMatrix<T> operator-(const T &scalar, const ygzMatrix<T> &matrix)
{
    T* result = new T[matrix.nElements];
    for (int i = 0; i < matrix.nElements; i++)
        result[i] = scalar - matrix.data[i];

    ygzMatrix<T> resultMatrix(matrix.nRows, matrix.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Matrix multiplication
template <class T>
ygzMatrix<T> operator*(const ygzMatrix<T> &lhs, const ygzMatrix<T> &rhs)
{
    if (lhs.nCols != rhs.nRows)
        throw std::invalid_argument("Number of columns in the first matrix must be equal to the number of rows in the second matrix");

    int nRows = lhs.nRows;
    int nCols = rhs.nCols;
    T* result = new T[nRows * nCols]; // Resulting matrix has dimensions (lhs.nRows, rhs.nCols)
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            int index = i * nCols + j;
            result[index] = 0; // Initialize element to 0
            for (int c = 0; c < lhs.nCols; c++) // for the ith row of the lhs
                for (int r = 0; r < rhs.nCols; r++) // for the jth column of the rhs
                    result[index] += lhs.getElement(i, c) * rhs.getElement(r, j); // Multiply and accumulate (dot product)
        }
    }

    ygzMatrix<T> resultMatrix(nRows, nCols, result);
    delete[] result;
    return resultMatrix;
}

// Matrix-scalar multiplication
template <class T>
ygzMatrix<T> operator*(const ygzMatrix<T> &matrix, const T &scalar)
{
    T* result = new T[matrix.nElements];
    for (int i = 0; i < matrix.nElements; i++)
        result[i] = matrix.data[i] * scalar;

    ygzMatrix<T> resultMatrix(matrix.nRows, matrix.nCols, result);
    delete[] result;
    return resultMatrix;
}

// Scalar-matrix multiplication
template <class T>
ygzMatrix<T> operator*(const T &scalar, const ygzMatrix<T> &matrix)
{
    return matrix * scalar; // Multiplication is commutative
}

/* ********************************************************************************************* */
/* 					    			  Checks                  									 */
/* ********************************************************************************************* */
template <class T>
bool ygzMatrix<T>::isSquare() const
{
    return nRows == nCols;
}

template <class T>
bool ygzMatrix<T>::closeEnough(T a, T b, double tolerance) const
{
    return fabs(a - b) < tolerance;
}

/* ********************************************************************************************* */
/* 					    			  Augmented Matrices    									 */
/* ********************************************************************************************* */
// Join two matrices into one (augmented matrix)
template <class T>
ygzMatrix<T> ygzMatrix<T>::join(const ygzMatrix<T> &m)
{
    if (nRows != m.nRows)
        throw std::invalid_argument("Matrices must have the same number of rows");

    T* newData = new T[nRows * (nCols + m.nCols)];
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < (nCols + m.nCols); j++)
        {
            int index = i * (nCols + m.nCols) + j;
            if (j < nCols)
                newData[index] = getElement(i, j);
            else
                newData[index] = m.getElement(i, j - nCols);
        }
    }
    ygzMatrix<T> result(nRows, (nCols + m.nCols), newData);
    delete[] newData;
    return result;
}

// Separate augmented matrix into two matrices
// Output is generated to m1 and m2 given as the input arguments
template <class T>
void ygzMatrix<T>::separate(ygzMatrix<T> *m1, ygzMatrix<T> *m2, int colNum) const
{
    if (colNum < 0 || colNum >= nCols)
        throw std::out_of_range("Column number out of range");

    m1->resize(nRows, colNum);
    m2->resize(nRows, nCols - colNum);

    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < colNum; j++)
            m1->setElement(i, j, getElement(i, j));
        for (int j = 0; j < nCols - colNum; j++)
            m2->setElement(i, j, getElement(i, colNum + j));
    }
}

/* ********************************************************************************************* */
/* 					    			  Gauss Jordan Elimination									 */
/* ********************************************************************************************* */
// Swaps row i with row j (in place)
template <class T>
void ygzMatrix<T>::swapRow(int i, int j)
{
    if (i == j)
        return;
    for (int c = 0; c < nCols; c++)
    {
        T temp = getElement(i, c);
        setElement(i, c, getElement(j, c));
        setElement(j, c, temp);
    }
}

// Multiplies row i by a factor (in place)
template <class T>
void ygzMatrix<T>::multRow(int i, T factor)
{
    for (int c = 0; c < nCols; c++)
        setElement(i, c, getElement(i, c) * factor);
}

// Multiplies row i by a factor and adds it to row j (in place)
template <class T>
void ygzMatrix<T>::multAddRow(int i, int j, T factor)
{
    for (int c = 0; c < nCols; c++)
        setElement(j, c, getElement(j, c) + factor * getElement(i, c));
}

// Function to find the row with the maximum value in a given column
// Returns the row index
template <class T>
int ygzMatrix<T>::findRowWithMax(int col, int startingRow) const
{
    T max = getElement(startingRow, col);
    int rowIndex = startingRow;
    for (int i = startingRow + 1; i < nRows; i++)
    {
        T val = getElement(i, col);
        if (val > max)
        {
            max = val;
            rowIndex = i;
        }
    }
    return rowIndex;
}

/* ********************************************************************************************* */
/* 					    			  Inverse Operations										 */
/* ********************************************************************************************* */
// Inverts the matrix in place
template <class T>
bool ygzMatrix<T>::inverseInPlace()
{
    if (!isSquare())
        throw std::invalid_argument("Matrix must be square to invert");

    // Create an identity matrix
    ygzMatrix<T> identityMatrix = ygzMatrix<T>::identity(nRows);
    ygzMatrix<T> augmentedMatrix = join(identityMatrix); // Augment the matrix with the identity matrix

    // Perform Gauss-Jordan elimination
    ygzMatrix<T>* lhs = new ygzMatrix<T>(*this);
    ygzMatrix<T>* res = new ygzMatrix<T>(nRows, nCols);

    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool complete = false;
    while ((!complete) && (count < maxCount))
    {
        for (int diagIndex = 0; diagIndex < nRows; diagIndex++)
        {
            cRow = diagIndex;
            cCol = diagIndex;
            int maxRow = augmentedMatrix.findRowWithMax(cCol, cRow); // Find row with max element in the cCol th column
            if (maxRow == -1)
                return false; // something went wrong
            augmentedMatrix.swapRow(cRow, maxRow); // Swap the row with the maximum value to the current row
            if (augmentedMatrix.getElement(cRow, cCol) != 1)
                augmentedMatrix.multRow(cRow, 1 / augmentedMatrix.getElement(cRow, cCol)); // Divide the row by the diagonal element to make it 1

            // Now process the rows and columns
            for (int r = cRow + 1; r < nRows; r++)
            {
                T currentElement = augmentedMatrix.getElement(r, cCol);
                T diagElement = augmentedMatrix.getElement(cRow, cCol);
                if (closeEnough(currentElement, 0.0) || closeEnough(diagElement, 0.0))
                    continue;

                T factor = -(currentElement / diagElement);
                augmentedMatrix.multAddRow(cRow, r, factor); // Make the elements below the diagonal element 0
            }
            for (int c = cCol + 1; c < nCols; c++)
            {
                T currentElement = augmentedMatrix.getElement(cRow, c);
                T diagElement = augmentedMatrix.getElement(cRow, cCol);
                if (closeEnough(currentElement, 0.0) || closeEnough(diagElement, 0.0))
                    continue;

                T factor = -(currentElement / diagElement);
                augmentedMatrix.multAddRow(c, cCol, factor); // Make the elements to the right of the diagonal element 0
            }
        }
        // Now seperate
        augmentedMatrix.separate(lhs, res, nCols);
        if (*lhs == identityMatrix) // we are done
        {
            complete = true;
            // res is now the inverse of the original matrix
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                    setElement(i, j, res->getElement(i, j));
            }
        }
        count++;
    }
    return complete;
}

template class ygzMatrix<int>;
template class ygzMatrix<long>;
template class ygzMatrix<float>;
template class ygzMatrix<double>;