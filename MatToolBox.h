/* EXAMPLE 10.3.1
   Copyright(c) 2005-10, S. D. Rajan
   Object-Oriented Numerical Analysis
   (c) James Kristoff 2011

  NOTES
   (1) Illustrates the development of a matrix toolbox
       based on the CVector and CMatrix classes
       developed earlier.
   (2) return value of most functions is either true or false.
       a true value signifies no error
       a false value indicates an input data error
*/
#ifndef JAMES_KRISTOFF_MATTOOLBOX_H
#define JAMES_KRISTOFF_MATTOOLBOX_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include "vectortemplate.h"
#include "matrixtemplate.h"

const int ELEMENTSPERLINE = 4;  // # of vector/matrix elements per line
const int FW = 10;              // field width

template <class T>
class CMatToolBox
{
    public:
        CMatToolBox ();
        ~CMatToolBox ();

        // vector-related functions
        void Display (const std::string& szMessage,
                      const CVector<T>& A) const;
        bool Add (const CVector<T>& A, const CVector<T>& B, 
                  CVector<T>& C);
        bool Subtract (const CVector<T>& A, const CVector<T>& B, 
                       CVector<T>& C);
        bool DotProduct (const CVector<T>& A,
                         const CVector<T>& B, T& product);
        bool Normalize (CVector<T>& A);
        void Scale (CVector<T>& A, T factor);
        T    MaxValue (const CVector<T>& A);
        T    MinValue (const CVector<T>& A);
        T    TwoNorm (const CVector<T>& A);
        T    MaxNorm (const CVector<T>& A);
        bool CrossProduct (const CVector<T>& A,
                           const CVector<T>& B, CVector<T>& C);

        // matrix-related functions
        void Display (const std::string& szMessage,
                      const CMatrix<T>& A) const;
        bool Add (const CMatrix<T>& A, const CMatrix<T>& B, 
                  CMatrix<T>& C);
        bool Subtract (const CMatrix<T>& A,
                       const CMatrix<T>& B, CMatrix<T>& C);
        bool Multiply (const CMatrix<T>& A,
                       const CMatrix<T>& B, CMatrix<T>& C);
        void Scale (CMatrix<T>& A, T scale);
        T    MaxNorm (const CMatrix<T>& A);
        bool Transpose (const CMatrix<T>& A,
                        CMatrix<T>& B);
        bool MatMultVec (const CMatrix<T>& A,
                         const CVector<T>& x,
                         CVector<T>& b);
        bool AxEqb (CMatrix<T>& A,
                    CVector<T>& x,
                    CVector<T>& b,
                    T TOL);
        bool LDLTFactorization (CMatrix<T>& A, T TOL);
        bool LDLTSolve (const CMatrix<T>& A, CVector<T>& x,
                        const CVector<T>& b);

        // helper function
        void GetFLOPStats (double& dAS, double& dM, double& dD) const;

    private:
        // these are declared double to avoid integer overflow
        double m_dASOP; // # of floating point additions and subtractions
        double m_dMOP;  // # of floating point multiplications
        double m_dDOP;  // # of floating point divisions

    protected:
};

// ctor
template <class T>
CMatToolBox<T>::CMatToolBox ()
// ==================================================================
// Function: default constructor
//    Input: none
//   Output: none
// ==================================================================
{
    m_dASOP = m_dMOP = m_dDOP = 0.0;
}

// dtor
template <class T>
CMatToolBox<T>::~CMatToolBox ()
// ==================================================================
// Function: destructor
//    Input: none
//   Output: none
// ==================================================================
{
}

// ---------------------------------------------------------------
// ------------------------ vector functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& szMessage,
                              const CVector<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: message string and the vector 
//   Output: None
// ==================================================================
{
    std::cout << '\n' << szMessage << '\n';
    std::cout.setf(std::ios::left);
    int n = A.GetSize();
    for (int i=1; i <= n; i++)
    {
        std::cout << "(" << i << ") "
                  << std::setw(FW) << A(i) << " ";
        if ((i % ELEMENTSPERLINE) == 0 || i == n)
            std::cout << '\n';
    }
}

template <class T>
bool CMatToolBox<T>::Add (const CVector<T>& A, const CVector<T>& B, 
                          CVector<T>& C)
// ==================================================================
// Function: adds two vectors and stores the result in the
//           third vector C = A + B
//    Input: vectors A, B, and C 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        return false;

    // add
    for (int i=1; i <= n; i++)
        C(i) = A(i) + B(i);
    m_dASOP += static_cast<double>(n);

    return true;
}

template <class T>
bool CMatToolBox<T>::Subtract (const CVector<T>& A, 
                               const CVector<T>& B, CVector<T>& C)
// ==================================================================
// Function: subtracts one vector from another and stores the result
//           in the third vector C = A - B
//    Input: vectors A, B, and C
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        return false;

    // subtract
    for (int i=1; i <= n; i++)
        C(i) = A(i) - B(i);
    m_dASOP += static_cast<double>(n);

    return true;
}

template <class T>
bool CMatToolBox<T>::DotProduct (const CVector<T>& A,
                                 const CVector<T>& B, T& product)
// ==================================================================
// Function: computes the dot product of two vectors such that
//           product = A dot B
//    Input: vectors A and B, variable to store the product
//   Output: product
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize())
        return false;

    // dot product
    product = 0;
    for (int i=1; i <= n; i++)
    {
        product += A(i) * B(i);
    }
    m_dASOP += static_cast<double>(n);
    m_dMOP += static_cast<double>(n);

    return true;
}

template <class T>
bool CMatToolBox<T>::Normalize (CVector<T>& A)
// ==================================================================
// Function: normalizes a vector
//    Input: vector A 
//   Output: normalized vector A 
// ==================================================================
{
    T sum = TwoNorm(A);
    int n = A.GetSize();

    if (sum == T(0))
        return false;
    else
    {
        for (int i=1; i <= n; i++)
            A(i) /= sum;
        m_dDOP += static_cast<double>(n);
    }

    return true;
}

template <class T>
void CMatToolBox<T>::Scale (CVector<T>& A, T c)
// ==================================================================
// Function: scales a vector by a constant c such that A = c A
//    Input: vector A and constant c 
//   Output: scaled vector A
// ==================================================================
{
    int n = A.GetSize();
    for (int i=1; i <= n; i++)
        A(i) *= c;

    m_dMOP += static_cast<double>(n);
}

template <class T>
T CMatToolBox<T>::MaxValue (const CVector<T>& A)
// ==================================================================
// Function: finds the largest value among all the elements in A
//    Input: vector A 
//   Output: return value is the largest element in A
// ==================================================================
{
    T max(0);
    int n = A.GetSize();

    for (int i=1; i <= n; i++)
    {
        if(A(i) > max)
            max = A(i);
    }

    return max;
}

template <class T>
T CMatToolBox<T>::MinValue (const CVector<T>& A)
// ==================================================================
// Function: finds the smallest value among all the elements in A
//    Input: vector A 
//   Output: return value is the smallest element in A
// ==================================================================
{
    T min(0);
    int n = A.GetSize();

    for (int i=1; i <= n; i++)
    {
        if(A(i) < min)
            min = A(i);
    }

    return min;
}

template <class T>
T CMatToolBox<T>::TwoNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two-norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{
    T sum(0);
    int n = A.GetSize();
    for (int i=1; i <= n; i++)
        sum += A(i) * A(i);

    m_dASOP += static_cast<double>(n);
    m_dMOP += static_cast<double>(n);

    return sqrt(sum);
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the inf-norm of vector A
//    Input: vector A 
//   Output: return value is the inf-norm
// ==================================================================
{
    T max(0);
    int n = A.GetSize();

    for (int i=1; i <= n; i++)
    {
        if(abs(A(i)) > max)
            max = abs(A(i));
    }

    return max;
}

template <class T>
bool CMatToolBox<T>::CrossProduct (const CVector<T>& A,
                                   const CVector<T>& B,
                                   CVector<T>& C)
// ==================================================================
// Function: computes the cross-product of two vectors and stores the
//           result in the third vector such that C = A x B
//           (3-dimensional space)
//    Input: vectors A, B and C
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != 3 || n != B.GetSize() || n != C.GetSize())
        return false;

    // compute the cross product
    C(1) = A(2) * B(3) - A(3) * B(2);
    C(2) = A(3) * B(1) - A(1) * B(3);
    C(3) = A(1) * B(2) - A(2) * B(1);

    m_dASOP += static_cast<double>(n);
    m_dMOP += static_cast<double>(2 * n);

    return true;
}

// ---------------------------------------------------------------
// ------------------------ matrix functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& szMessage,
                              const CMatrix<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a matrix
//           rowwise
//    Input: message string and the matrix
//   Output: None
// ==================================================================
{
    std::cout << '\n' << szMessage << '\n';
    std::cout.setf(std::ios::left);
    int R = A.GetRows();
    int C = A.GetColumns();

    for (int i=1; i <= R; i++)
    {
        for (int j=1; j <= C; j++)
        {
            std::cout << "(" << i << "," << j << ") "
                      << std::setw(FW) << A(i,j) << " ";
            if ((j % ELEMENTSPERLINE) == 0 || (i == R && j == C))
                std::cout << '\n';
        }
        std::cout << '\n';
    }
}

template <class T>
bool CMatToolBox<T>::Add (const CMatrix<T>& A, const CMatrix<T>& B, 
                          CMatrix<T>& C)
// ==================================================================
// Function: adds two matrices and stores the result in the
//           third matrix C = A + B
//    Input: matrices A, B, and C 
//   Output: matrix C
// ==================================================================
{
    int nR = A.GetRows();
    int nC = A.GetColumns();

    //Check for incompatible matrices
    if (nR != B.GetRows() || nR != C.GetRows() ||
        nC != B.GetColumns() || nC != C.GetColumns())
        return false;

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            C(i,j) = A(i,j) + B(i,j);
        }
    }

    m_dASOP += static_cast<double>(nR * nC);

    return true;
}

template <class T>
bool CMatToolBox<T>::Subtract (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: subtracts one matrix from another and stores the result
//           in the third matrix C = A - B
//    Input: matrices A, B, and C
//   Output: matrix C
// ==================================================================
{
    int nR = A.GetRows();
    int nC = A.GetColumns();

    //Check for incompatible matrices
    if (nR != B.GetRows() || nR != C.GetRows() ||
        nC != B.GetColumns() || nC != C.GetColumns())
        return false;

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            C(i,j) = A(i,j) - B(i,j);
        }
    }

    m_dASOP += static_cast<double>(nR * nC);

    return true;
}

template <class T>
bool CMatToolBox<T>::Multiply (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: multiplies two matrices and stores the result
//           in the third matrix C = A * B
//    Input: matrices A, B, and C
//   Output: matrix C
// ==================================================================
{
    int nRa = A.GetRows();
    int nCa = A.GetColumns();
    int nCb = B.GetColumns();

    //Check for incompatible matrices
    if (nCa != B.GetRows() || nRa != C.GetRows() || 
        nCb != C.GetColumns())
        return false;

    C.Set(0);

    for (int j=1; j <= nCb; j++)
    {
        for (int i=1; i <= nRa; i++)
        {
            for (int k=1; k <= nCa; k++)
            {
                C(i,j) += A(i,k) * B(k,j);
            }
        }
    }

    m_dASOP += static_cast<double>(nRa * nCa * nCb);
    m_dMOP += static_cast<double>(nRa * nCa * nCb);

    return true;
}

template <class T>
void CMatToolBox<T>::Scale (CMatrix<T>& A, T c)
// ==================================================================
// Function: scales all the elements of a matrix by a constant c
//           such that A = c A
//    Input: matrix A and constant c
//   Output: scaled matrix A
// ==================================================================
{
    int nR = A.GetRows();
    int nC = A.GetColumns();

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            A(i,j) *= c;
        }
    }

    m_dMOP += static_cast<double>(nR * nC);
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CMatrix<T>& A)
// ==================================================================
// Function: computes the max norm of matrix A
//    Input: matrix A 
//   Output: return value is the max norm
// ==================================================================
{
    T max(0);
    int nR = A.GetRows();
    int nC = A.GetColumns();

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            if(abs(A(i,j)) > max)
                max = abs(A(i,j));
        }
    }

    return max;
}

template <class T>
bool CMatToolBox<T>::Transpose (const CMatrix<T>& A, CMatrix<T>& B)
// ==================================================================
// Function: computes the transpose of a matrix and stores the result
//           in another matrix B = A(T)
//    Input: matrices A and B
//   Output: matrix B
// ==================================================================
{
    int nR = A.GetRows();
    int nC = A.GetColumns();

    if (nR != B.GetColumns() || nC != B.GetRows())
        return false;

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            B(j,i) = A(i,j);
        }
    }

    return true;
}

template <class T>
bool CMatToolBox<T>::MatMultVec (const CMatrix<T>& A,
                                 const CVector<T>& x,
                                 CVector<T>& b)
// ==================================================================
// Function: multiplies a matrix and a vector and stores the result
//           in a vector b = A * x
//    Input: matrix A and vectors x and b
//   Output: vector b
// ==================================================================
{
    int nR = A.GetRows();
    int nC = A.GetColumns();

    if (nC != x.GetSize() || nR != b.GetSize())
        return false;

    b.Set(0);

    for (int i=1; i <= nR; i++)
    {
        for (int j=1; j <= nC; j++)
        {
            b(i) += A(i,j) * x(j);
        }
    }

    m_dASOP += static_cast<double>(nR * nC);
    m_dMOP += static_cast<double>(nR * nC);

    return true;
}

template <class T>
bool CMatToolBox<T>::AxEqb (CMatrix<T>& A,
                            CVector<T>& x,
                            CVector<T>& b,
                            T TOL)
// ==================================================================
// Function: solves A x = b using Gaussian Elimination Technique
//           (this version only for one rhs vector)
//    Input: Matrix A, Vectors x and b, and a constant tolerance
//   Output: Vector x
//           A false return value indicates a singular A matrix.
// ==================================================================
{
    int i, j, k;   // loop indices
    T c;           // multiplier (Step 4)

    // number of equations to solve
    int n = A.GetRows();

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)              // Step 1
    {
        if (fabs(A(k,k)) <= TOL)          // Step 2
            return false;
        for (i=k+1; i <= n; i++)          // Step 3
        {
            c = A(i,k)/A(k,k);            // Step 4
            for (j=k+1; j <= n; j++)      // Step 5
                A(i,j) -= c * A(k,j);     // Step 6
            x(i) -= c * x(k);             // Step 8
        }                                 // Step 9
        int nC = n-k;
        if (nC > 0)
        {
            m_dDOP += static_cast<double>(nC);
            m_dMOP += static_cast<double>(nC*nC);
            m_dASOP += static_cast<double>(nC+nC*nC);
        }
    }                                     // Step 10 

    // back substitution
    if (fabs(A(n,n)) <= TOL)              
        return false;
    x(n) /= A(n,n);                       // Step 11

    for (i=n-1; i >= 1; i--)              // Step 12
    {
        T sum = T(0);
        for (j=i+1; j <= n; j++)
            sum += A(i,j) * x(j);         // Step 13
        if ((n-i) > 0)
        {
            m_dASOP += static_cast<double>(n-i);
            m_dMOP += static_cast<double>(n-i);
        }
        x(i) = (x(i) - sum)/A(i,i);       // Step 14
    }                                     // Step 15
    m_dASOP += static_cast<double>(n);
    m_dDOP += static_cast<double>(n+1);

    return true;
}

template <class T>
bool CMatToolBox<T>::LDLTFactorization (CMatrix<T>& A,
                                        T TOL)
// ==================================================================
// Function: carries out LDL(T) factorization of matrix A
//           A is replaced with L and D. A is a symmetric matrix.
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
    // number of equations to solve
    int n = A.GetRows();

    for(int i = 1; i <= n; i++)                     // Step 1
    {
        for(int j = 1; j <= i-1; j++)
        {
            A(i,i) -= A(i,j) * A(i,j) * A(j,j);     // Step 2
            if(A(i,i) < TOL)
            {
                m_dASOP += static_cast<double>(j);
                m_dMOP += static_cast<double>(2.0 * j);
                return false;
            }
        }
        m_dASOP += static_cast<double>(i-1);
        m_dMOP += static_cast<double>(2.0 * (i - 1.0));
        for(int j = i+1; j <= n; j++)
        {
            for(int k = 1; k <= i-1; k++)
                A(j,i) -= A(j,k) * A(k,k) * A(i,k);
            A(j,i) /= A(i,i);                       // Step 3
            m_dASOP += static_cast<double>(i-1);
            m_dMOP += static_cast<double>(2.0 * (i - 1.0));
            m_dDOP += 1.0;
        }
    }                                               // Step 4

    return true;
}

template <class T>
bool CMatToolBox<T>::LDLTSolve (const CMatrix<T>& A,
                                CVector<T>& x,
                                const CVector<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and D terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    // number of equations to solve
    int n = A.GetRows();

    x = b;
    for(int i = 2; i <= n; i++)
    {
        for(int j = 1; j <= i - 1; j++)
        {
            x(i) -= A(i,j) * x(j);                  // Step 6
        }
        m_dASOP += static_cast<double>(i-1);
        m_dMOP += static_cast<double>(i-1);
    }
    for(int i = n; i >= 1; i--)
    {
        x(i) /= A(i,i);                         // Step 7
        for(int j = i+1; j <= n; j++)
        {
            x(i) -= A(j,i) * x(j);                  // Step 8
        }
        m_dASOP += static_cast<double>(n - i);
        m_dMOP += static_cast<double>(n - i);
        m_dDOP += 1.0;
    }
    return true;
}

template <class T>
void CMatToolBox<T>::GetFLOPStats (double& dAS, double& dM, 
                                   double& dD) const
// ==================================================================
// Function: retrieves floating point operations
//    Input: variables to store +-, * and / operations
//   Output: variables with their values
// ==================================================================
{
    dAS = m_dASOP;
    dM  = m_dMOP;
    dD  = m_dDOP;
}

#endif
