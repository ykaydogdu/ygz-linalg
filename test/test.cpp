#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector> // Add this line to include the vector header

using namespace std;

#include "ygzlinalg.hpp"

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

int main()
{
    cout.precision(10);

    cout << "Test starting\n";

    fstream fin;
    fin.open("./test/test_inverse.csv", ios::in);
    
    // first line contains the number of rows followed by '# '
    string line, word, temp;
    getline(fin, line);
    stringstream ss(line);
    getline(ss, word, ' ');
    getline(ss, word, ' ');
    int n = stoi(word);
    if (n <= 0)
    {
        cout << "Invalid number of rows\n";
        return 1;
    }

    // read the matrix
    ygzMatrix<double>* X = new ygzMatrix<double>[2 * n];
    ygzMatrix<double>* A = new ygzMatrix<double>[n];
    ygzMatrix<double>* A_inv = new ygzMatrix<double>[n];
    for (int i = 0; i < 2 * n; i++)
    {
        X[i].resize(n, n);
    }

    // read csv
    int i = 0;
    while (getline(fin, line))
    {
        // read each element of the row
        stringstream ss(line);
        int j = 0;
        while (getline(ss, word, ','))
        {
            if (j < n)
            {
                X[((int)i/n)].setElement(i%n, j, stod(word));
            }
            else
            {
                X[n + ((int)i/n)].setElement(i%n, j-n, stod(word));
            }
            j++;
        }
        i++;
    }

    // close the file
    fin.close();

    // split the matrix into A and A_inv
    for (int i = 0; i < 2 * n; i++)
    {
        if (i < n)
        {
            A[i] = X[i];
        }
        else
        {
            A_inv[i - n] = X[i];
        }
    }

    // test the inverse function
    int passedCount = 0;
    int failedCount = 0;
    for (int i = 0; i < n; i++)
    {
        bool success = false;
        A[i].inverseInPlace();
        if (A[i] == A_inv[i])
        {
            passedCount++;
            success = true;
        }
        else
        {
            failedCount++;
        }
        cout << "Test " << i << ": " << success << endl;
    }

    cout << "Passed: " << passedCount << endl;
    cout << "Failed: " << failedCount << endl;
    return 0;
}