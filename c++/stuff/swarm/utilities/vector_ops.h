/* Defines a set of standard mathematical "operators" on mathematical
 * "vectors"
 */
#ifndef  __VECTOR_OPS_H__
#define __VECTOR_OPS_H__

#include <vector>
#include <functional>
#include <numeric>
#include <algorithm>
#include <iostream>

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
T dotProd(const std::vector<T>& v1, std::vector<T>& v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), T(0.0));
}

template <typename T>
T norm(const std::vector<T>& v1) {
    return sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), T(0.0)));
}

#endif
