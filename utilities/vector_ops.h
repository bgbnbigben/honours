/* Defines a set of standard mathematical "operators" on mathematical
 * "vectors"
 */
#ifndef  __VECTOR_OPS_H__
#define __VECTOR_OPS_H__

#include <functional>
#include <numeric>
#include <algorithm>
#include <iostream>

#include <vector>

/* Addition operators */
template <typename T>
std::vector<T> operator+(const std::vector<T> &v1, const std::vector<T> &v2) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), ret.begin(), std::plus<T>());
    return ret;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &v1, const T& c1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind2nd(std::plus<T>(), c1));
    return ret;
}

template <typename T>
std::vector<T> operator+(const T& c1, const std::vector<T> &v1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind1st(std::plus<T>(), c1));
    return ret;
}

template <typename T>
void operator+=(std::vector<T> &v1, const std::vector<T> &v2) {
    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::plus<T>());
}

template <typename T>
void operator+=(std::vector<T> &v1, const T &c1) {
    std::transform(v1.begin(), v1.end(), v1.begin(), std::bind2nd(std::plus<T>(), c1));
}


/* Subtraction operators and the negation operator */
template <typename T>
std::vector<T> operator-(std::vector<T> &v1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::negate<T>());
    return ret;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &v1, const std::vector<T> &v2) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), ret.begin(), std::minus<T>());
    return ret;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &v1, const T& c1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind2nd(std::minus<T>(), c1));
    return ret;
}

template <typename T>
std::vector<T> operator-(const T& c1, const std::vector<T> &v1) {
    return -v1 + c1; // TODO: Check optimisations.
}

template <typename T>
void operator-=(std::vector<T> &v1, const std::vector<T> &v2) {
    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::minus<T>());
}

template <typename T>
void operator-=(std::vector<T> &v1, const T &c1) {
    std::transform(v1.begin(), v1.end(), v1.begin(), std::bind2nd(std::minus<T>(), c1));
}

/* Multiplication operators */
/* Scales a vector by another vector element-wise */
template <typename T>
std::vector<T> operator*(const std::vector<T> &v1, const std::vector<T> &v2) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), ret.begin(), std::multiplies<T>());
    return ret;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &v1, const T& c1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind2nd(std::multiplies<T>(), c1));
    return ret;
}

template <typename T>
std::vector<T> operator*(const T& c1, const std::vector<T> &v1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind1st(std::multiplies<T>(), c1));
    return ret;
}

template <typename T>
void operator*=(std::vector<T> &v1, const std::vector<T> &v2) {
    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::multiplies<T>());
}

template <typename T>
void operator*=(std::vector<T> &v1, const T &c1) {
    std::transform(v1.begin(), v1.end(), v1.begin(), std::bind2nd(std::multiplies<T>(), c1));
}

/* Division operators */

/* Scales a vector by another vector element-wise */
template <typename T>
std::vector<T> operator/(const std::vector<T> &v1, const std::vector<T> &v2) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), ret.begin(), std::divides<T>());
    return ret;
}

template <typename T>
std::vector<T> operator/(const std::vector<T> &v1, const T& c1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind2nd(std::divides<T>(), c1));
    return ret;
}

template <typename T>
std::vector<T> operator/(const T& c1, const std::vector<T> &v1) {
    std::vector<T> ret(v1.size());
    std::transform(v1.begin(), v1.end(), ret.begin(), std::bind1st(std::divides<T>(), c1));
    return ret;
}

template <typename T>
void operator/=(std::vector<T> &v1, const std::vector<T> &v2) {
    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::divides<T>());
}

template <typename T>
void operator/=(std::vector<T> &v1, const T &c1) {
    std::transform(v1.begin(), v1.end(), v1.begin(), std::bind2nd(std::divides<T>(), c1));
}

/**
 * Overload the output stream function.
 */
template <typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &v) {
	int N = v.size();
	out << "size: " << N << " x 1" << "\n";

	for( int i=0; i<N; ++i )
		out << v[i] << "\n";

	return out;
}

template <typename T>
T dotProd(const std::vector<T> v1, std::vector<T> v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), T(0.0));
}

template <typename T>
T norm(const std::vector<T>& v1) {
    return sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), T(0.0)));
}

namespace std {
    template class vector<double>;
}

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

template <typename T>
std::vector<T> inverse(std::vector<T> mat) {
    int N = mat.size();
    assert(N && N == mat[0].size());
    double* A = new double[N*N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i*(j+1)] = mat[i][j];
        }
    }
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mat[i][j] = A[i*(j+1)];
        }
    }
    return mat;
}

template <typename T>
std::vector<T> transpose(std::vector<T> mat) {
    int N = mat.size();
    assert(N && N == mat[0].size());
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            auto ij = mat[i][j];
            auto ji = mat[j][i];
            mat[i][j] = ji;
            mat[j][i] = ij;
        }
    }

    return mat;
}

/* matrix mult */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& m, const std::vector<T> &v) {
    assert(m.size() == v.size());
    std::vector<T> ret(v.size());
    for (int i = 0; i < v.size(); i++) {
        auto temp = m[i];
        ret[i] = dotProd<double>(temp, v);
    }
    return ret;
}



#endif
