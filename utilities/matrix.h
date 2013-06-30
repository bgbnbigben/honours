#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <cassert>
#include "vector_ops.h"

template <typename T>
class Matrix {

    public:
        // constructors and destructor
        Matrix() : pv0(0), prow0(0), nRow(0), nColumn(0), nTotal(0) {};
        Matrix(const Matrix<T> &A);
        Matrix(int rows, int columns, const T &x = T(0));
        Matrix(int rows, int columns, const T *v);
        ~Matrix();

        // assignments
        Matrix<T>& operator=(const Matrix<T> &A);
        Matrix<T>& operator=(const T &x);

        // accessors
        T* operator[](int i);
        const T* operator[](int i) const;
        T& operator()(int row, int column);
        const T& operator()(int row, int column) const;

        // type conversion
        operator T*();
        operator const T*() const;

        // others
        long size() const;
        int dim(int dimension) const;
        int rows() const;
        int cols() const;
        Matrix<T>& resize(int rows, int columns);
        std::vector<T> getRow(int row) const;
        std::vector<T> getColumn(int column) const;
        void setRow(const std::vector<T> &v, int row);
        void setColumn(const std::vector<T> &v, int column);

        // computed assignment
        Matrix<T>& operator+=(const T&);
        Matrix<T>& operator+=(const Matrix<T>&);
        Matrix<T>& operator-=(const T&);
        Matrix<T>& operator-=(const Matrix<T>&);
        Matrix<T>& operator*=(const T&);
        // WARNING: element-by-element
        Matrix<T>& operator*=(const Matrix<T>&);
        Matrix<T>& operator/=(const T&);
        // WARNING: element-by-element
        Matrix<T>& operator/=(const Matrix<T>&);

    private:

        // Full data
        T *pv0;

        // row pointer's pointer
        T **prow0;

        // row number, column number and total number
        int	 nRow;
        int	 nColumn;
        long nTotal;

        void init(int rows, int columns);
        void copyFromArray(const T *v);
        void setByScalar(const T &x);
        void destroy();

};
// class Matrix

// input and output
template<typename T>
std::ostream& operator<<( std::ostream&, const Matrix<T>& );
template<typename T>
std::istream& operator>>( std::istream&, Matrix<T>& );

// arithmetic operators
template<typename T>
Matrix<T> operator- (const Matrix<T>&);
template<typename T>
Matrix<T> operator+ (const Matrix<T>&, const T&);
template<typename T>
Matrix<T> operator+ (const T&, const Matrix<T>&);
template<typename T>
Matrix<T> operator+ (const Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T> operator- (const Matrix<T>&, const T&);
template<typename T>
Matrix<T> operator- (const T&, const Matrix<T>&);
template<typename T>
Matrix<T> operator- (const Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T> operator* (const Matrix<T>&, const T&);
template<typename T>
Matrix<T> operator* (const T&, const Matrix<T>&);
template<typename T>
Matrix<T> operator* (const Matrix<T>&, const Matrix<T>&);
template<typename T>
std::vector<T> operator* (const Matrix<T>&, const std::vector<T>&);
template<typename T>
Matrix<T> operator/ (const Matrix<T>&, const T&);
template<typename T>
Matrix<T> operator/ (const T&, const Matrix<T>&);

// optimizied version of multiplication
template<typename T>
Matrix<T>& optMult (const Matrix<T>&, const Matrix<T>&, Matrix<T>&);
template<typename T>
std::vector<T>& optMult (const Matrix<T>&, const std::vector<T>&, std::vector<T>&);

// element-by-element multiplication and division
template<typename T>
Matrix<T> elemMult (const Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T> elemDivd (const Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T>& elemMultEq (Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T>& elemDivdEq (Matrix<T>&, const Matrix<T>&);

// transpose and conjugate transpose
template<typename T> Matrix<T> trT (const Matrix<T>&);
template<typename T> Matrix<T> trH (const Matrix<T>&);

// transpose multiplication
template<typename T>
Matrix<T> trMult (const Matrix<T>&, const Matrix<T>&);
template<typename T>
std::vector<T> trMult (const Matrix<T>&, const std::vector<T>&);
template<typename T>
Matrix<T> multTr (const Matrix<T>&, const Matrix<T>&);
template<typename T>
Matrix<T> multTr (const std::vector<T>&, const std::vector<T>&);

template<typename T> Matrix<T> eye (int, const T&);

// utilities
template<typename T> T norm (const Matrix<T>&);

/**
 * initialize
 */
template <typename T>
void Matrix<T>::init (int rows, int columns) {
	nRow = rows;
	nColumn = columns;
	nTotal = nRow * nColumn;

	pv0 = new T[nTotal];
	prow0 = new T*[nRow];

	T *p = pv0;
	for (int i = 0; i < nRow; i++) {
		prow0[i] = p;
		p += nColumn;
	}
}

/**
 * copy matrix from normal array
 */
template <typename T>
inline void Matrix<T>::copyFromArray( const T *v ) {
    std::copy(v, v+nTotal, pv0);
}

/**
 * set matrix by a scalar
 */
template <typename T>
inline void Matrix<T>::setByScalar( const T &x ) {
    std::fill(pv0, pv0+nTotal, x);
}

/**
 * destroy the matrix
 */
template <typename T>
void Matrix<T>::destroy() {
	if (pv0 == NULL)
		return;
	else
		delete []pv0;

	if (prow0 != NULL)
		delete []prow0;
}

/**
 * constructors and destructor
 */
template <typename T>
Matrix<T>::Matrix (const Matrix<T> &A) {
	init(A.nRow, A.nColumn);
	copyFromArray(A.pv0);
}

template <typename T>
Matrix<T>::Matrix (int rows, int columns, const T &x) {
	init(rows,columns);
	setByScalar(x);
}

template <typename T>
Matrix<T>::Matrix(int rows, int columns, const T *arrays) {
	init(rows, columns);
	copyFromArray(arrays);
}

template <typename T>
Matrix<T>::~Matrix() {
	destroy();
}

/**
 * overload evaluate operator = from matrix to matrix
 */
template <typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T> &A) {
	if (pv0 == A.pv0)
		return *this;

	if (nRow == A.nRow && nColumn == A.nColumn)
		copyFromArray (A.pv0);
	else {
		destroy();
		init (A.nRow, A.nColumn);
		copyFromArray (A.pv0);
	}

	return *this;
}

/**
 * overload evaluate operator = from scalar to matrix
 */
template <typename T>
inline Matrix<T>& Matrix<T>::operator= (const T &x) {
	setByScalar(x);
	return *this;
}


/**
 * overload operator [] for 0-offset access
 */
template <typename T>
inline T* Matrix<T>::operator[] (int i) {
	return prow0[i];
}

template <typename T>
inline const T* Matrix<T>::operator[] (int i) const {
	return prow0[i];
}

/* 1-offset access */
template <typename T>
inline T& Matrix<T>::operator()(int row, int column) {
    return prow0[row-1][column-1];
}

template <typename T>
inline const T& Matrix<T>::operator() (int row, int column) const {
    return prow0[row-1][column-1];
}


/**
 * type conversion functions
 */
template <typename T>
inline Matrix<T>::operator T*() {
	return pv0;
}

template <typename T>
inline Matrix<T>::operator const T*() const {
	return pv0;
}

/**
 * get the matrix's size
 */
template <typename T>
inline long Matrix<T>::size() const {
	return nTotal;
}

template <typename T>
inline int Matrix<T>::rows() const {
    return nRow;
}

template <typename T>
inline int Matrix<T>::cols() const {
    return nColumn;
}
/**
 * compound assignment operators +=
 */
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const T &x) {
    std::transform(pv0, pv0+nTotal, pv0, std::bind2nd(std::plus<T>(), x));
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &rhs) {
    assert( nRow == rhs.rows() );
    assert( nColumn == rhs.cols() );
    std::transform(pv0, pv0+nTotal, rhs.pv0, pv0, std::plus<T>());
    return *this;
}


template<typename T>
inline Matrix<T> operator+(const Matrix<T> &A, const T &x) {
	Matrix<T> tmp( A );
	return tmp += x;
}

template<typename T>
inline Matrix<T> operator+(const T &x, const Matrix<T> &A) {
	return A + x;
}

template<typename T>
inline Matrix<T> operator+(const Matrix<T> &A1, const Matrix<T> &A2) {
	Matrix<T> tmp(A1);
	return tmp += A2;
}

/* Subtraction */

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const T &x) {
    std::transform(pv0, pv0+nTotal, pv0, std::bind2nd(std::minus<T>(), x));
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &rhs) {
    assert(nRow == rhs.rows());
    assert(nColumn == rhs.cols());
    std::transform(pv0, pv0+nTotal, rhs.pv0, pv0, std::minus<T>());
    return *this;
}

template<typename T>
Matrix<T> operator-(const Matrix<T> &A) {
    int rows = A.rows(), columns = A.cols();

	Matrix<T> tmp(rows, columns);
    // TODO: lol
	for(int i = 0; i < rows; ++i)
		for(int j = 0; j < columns; ++j)
			tmp[i][j] = -A[i][j];

	return tmp;
}

template<typename T>
inline Matrix<T> operator-(const Matrix<T> &A, const T &x) {
	Matrix<T> tmp(A);
	return tmp -= x;
}

template<typename T>
inline Matrix<T> operator-(const T &x, const Matrix<T> &A) {
	Matrix<T> tmp(A);
	return -tmp += x;
}

template<typename T>
inline Matrix<T> operator-(const Matrix<T> &A1, const Matrix<T> &A2) {
	Matrix<T> tmp(A1);
	return tmp -= A2;
}

/* Multiplication */
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T &x) {
    std::transform(pv0, pv0+nTotal, pv0, std::bind2nd(std::multiplies<T>(), x));
    return *this;
}

// WARNING: this is element-by-element multiplication
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &rhs) {
    assert(nRow == rhs.rows());
    assert(nColumn == rhs.cols());
    std::transform(pv0, pv0+nTotal, rhs.pv0, pv0, std::multiplies<T>());
    return *this;
}

template <typename T>
inline Matrix<T> operator*(const Matrix<T> &A, const T &x) {
	Matrix<T> tmp(A);
	return tmp *= x;
}

template <typename T>
inline Matrix<T> operator*(const T &x, const Matrix<T> &A) {
	return A * x;
}

template <typename T>
Matrix<T> operator*(const Matrix<T> &A1, const Matrix<T> &A2) {
	assert(A1.cols() == A2.rows());

	int rows = A1.rows();
	int columns = A2.cols();

	Matrix<T> tmp(rows, columns);
    mult(A1, A2, tmp);

	return tmp;
}


template <typename T>
std::vector<T> operator*( const Matrix<T> &A, const std::vector<T> &b) {
	assert(A.cols() == b.size());

	int rows = A.rows();

	std::vector<T> tmp(rows);
    mult(A, b, tmp);

	return tmp;
}


/* Division */
template <typename T>
Matrix<T>& Matrix<T>::operator/=(const T &x) {
    std::transform(pv0, pv0+nTotal, pv0, std::bind2nd(std::divides<T>(), x));
    return *this;
}

// WARNING: this is element-by-element division
template <typename T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T> &rhs) {
    assert(nRow == rhs.rows());
    assert(nColumn == rhs.cols());
    std::transform(pv0, pv0+nTotal, rhs.pv0, pv0, std::divides<T>());
    return *this;
}

template <typename T>
inline Matrix<T> operator/(const Matrix<T> &A, const T &x) {
	Matrix<T> tmp(A);
	return tmp /= x;
}

template <typename T>
Matrix<T> operator/(const T &x, const Matrix<T> &A) {
	int rows = A.rows();
	int clumns = A.cols();

	Matrix<T> tmp(rows,clumns);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < clumns; ++j)
			tmp[i][j] = x / A[i][j];

	return tmp;
}

template <typename T>
std::ostream& operator<<(std::ostream &out, const Matrix<T> &A) {
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < columns; ++j)
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}

/**
 * This is an optimized version of matrix multiplication,
 * where the destination matrix has already been allocated.
 */
template <typename T>
Matrix<T>& mult(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C) {
    int M = A.rows();
    int N = B.cols();
    int K = A.cols();

    assert( B.rows() == K );

    C.resize(M, N);
    T        sum;
    const T  *pRow,
             *pCol;

    for (int i = 0; i < M; i++) 
        for (int j = 0; j < N; ++j) {
            pRow  = &A[i][0];
            pCol  = &B[0][j];
            sum = 0;

            for (int k = 0; k < K; ++k) {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol += N;
            }
            C[i][j] = sum;
        }
    return C;
}

template <typename T>
std::vector<T>& mult(const Matrix<T> &A, const std::vector<T> &b, std::vector<T> &c) {
    int M = A.rows();
    int N = A.cols();

    assert( b.size() == N );

    c.resize(M);
    T        sum;
    const T  *pRow,
             *pCol;

    for (int i = 0; i < M; i++)  {
        pRow  = &A[i][0];
        pCol  = &b[0];
        sum = 0;

        for (int j = 0; j < N; ++j)  {
            sum += (*pRow) * (*pCol);
            pRow++;
            pCol++;
        }
        c[i] = sum;
    }
    return c;
}

/**
 * matrix-matrix elementwise operations
 */
template<typename T>
inline Matrix<T> elemMult (const Matrix<T> &A1, const Matrix<T> &A2) {
	Matrix<T> tmp (A1);
	return tmp *= A2;
}

template <typename T>
inline Matrix<T>& elemMultEq (Matrix<T> &A1, const Matrix<T> &A2) {
    return A1 *= A2;
}

template <typename T>
inline Matrix<T> elemDivd (const Matrix<T> &A1, const Matrix<T> &A2) {
	Matrix<T> tmp (A1);
	return tmp /= A2;
}

template <typename T>
inline Matrix<T>& elemDivdEq(Matrix<T> &A1, const Matrix<T> &A2) {
    return A1 /= A2;
}

/**
 * matrix tranpose
 */
template <typename T>
Matrix<T> trT(const Matrix<T> &A) {
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<T> tmp(rows, clumns);
	for (int i = 0; i < rows; ++i) 
		for (int j = 0; j < clumns; ++j) 
			tmp[i][j] = A[j][i];

	return tmp;
}

/**
 * matrix conjugate tranpose
 */
template <typename T>
Matrix<T> trH(const Matrix<T> &A) {
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<T> tmp(rows, clumns);
	for (int i = 0; i < rows; ++i) 
		for (int j = 0; j < clumns; ++j) 
			tmp[i][j] = conj(A[j][i]);

	return tmp;
}

/**
 * matrix-matrix tranpose multiplication: A^T * B.
 */
template <typename T>
Matrix<T> trMult(const Matrix<T> &A1, const Matrix<T> &A2) {
	assert( A1.rows() == A2.rows() );

	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<T> tmp(rows, columns);
    for (int i = 0; i < rows; ++i) 
		for (int j = 0; j < columns; ++j) 
			for (int k = 0; k < K; ++k) 
			   tmp[i][j] += A1[k][i] * A2[k][j];

	return tmp;
}

/**
 * matrix-vector tranpose multiplication: A^T * b.
 */
template <typename T>
std::vector<T> trMult(const Matrix<T> &A, const std::vector<T> &v) {
	assert( A.rows() == v.size() );

	int rows = A.rows();
	int columns = A.cols();

	std::vector<T> tmp(columns);
    for (int i = 0; i < columns; ++i) 
		for (int j = 0; j < rows; ++j) 
			tmp[i] +=  A[j][i] * v[j];

	return tmp;
}

/**
 * matrix-matrix tranpose multiplication: A * B^T.
 */
template <typename T>
Matrix<T> multTr(const Matrix<T> &A1, const Matrix<T> &A2) {
	assert( A1.cols() == A2.cols() );

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<T> tmp(rows, columns);
    for (int i = 0; i < rows; ++i) 
		for (int j = 0; j < columns; ++j) 
			for (int k = 0; k < K; ++k) 
			   tmp[i][j] += A1[i][k] * A2[j][k];

	return tmp;
}

/**
 * vector-vector tranpose multiplication: a * b^T.
 */
template <typename T>
Matrix<T> multTr(const std::vector<T> &a, const std::vector<T> &b) {
	int rows = a.size();
	int columns = b.size();

	Matrix<T> tmp(rows, columns);
	for (int i = 0; i < rows; ++i) 
		for (int j = 0; j < columns; ++j) 
			tmp[i][j] = a[i]*b[j];

	return tmp;
}

/**
 * Generate the identity matrix.
 */
template <typename T>
Matrix<T> eye(int N, const T &x) {
    Matrix<T> tmp(N, N);
	for (int i = 0; i < N; ++i) 
		tmp[i][i] = x;

	return tmp;
}

/**
 * Compute Frobenius norm of matrix.
 */
template <typename T>
T norm(const Matrix<T> &A) {
	int m = A.rows();
	int n = A.cols();

	T sum = 0;
	for (int i = 1; i <= m; ++i) 
		for (int j = 1; j <= n; ++j) 
            sum += A(i,j) * A(i,j);

	return sqrt(sum);
}
#endif
