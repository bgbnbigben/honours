#ifndef LINESEARCH_H
#define LINESEARCH_H

#include <utilities/vector_ops.h>
#include <vector>
#include <algorithm>

template <typename T, typename F>
class LineSearch {
    public:

        LineSearch() : funcNum_(0), success_(true) {};
        ~LineSearch() {};

        T getStep(F &f, std::vector<T> &xk, std::vector<T> &dk, int maxItr=10);

        int numFuncEvaluations() const;
        bool isSuccess() const;

    protected:
        int     funcNum_;
        bool    success_;

};

/* x_k = point, d_k = direction, function = f */
template <typename T, typename F>
T LineSearch<T, F>::getStep(F &f, std::vector<T> &xk, std::vector<T> &dk, int maxItr) {
    // Set line search parameters that everyone uses.
    T mu = T(0.001),
          kUp = T(0.5),
          kLow = T(0.1),
          alpha = T(1.0),
          alphaMin,
          alphaMax;

    T fNew, fk = f(xk);

    std::vector<T> xNew, gk = f.grad(xk);

    T gd = dotProd(gk, dk);

    for (int i = 0; i < maxItr; ++i) {
        xNew = xk + alpha*dk;
        fNew = f(xNew);
        this->funcNum_++;

        if (fNew < fk + mu*alpha*gd) {
            this->success_ = true;
            return alpha;
        } else {
            alphaMin = kLow*alpha;
            alphaMax = kUp*alpha;

            // Compute the step by using quadratic polynomial interpolation.
            alpha = T(-0.5) * alpha * alpha * gd / (fNew - fk - alpha * gd);
            alpha = std::max(alphaMin, std::min(alphaMax, alpha));
        }
    }

    if (fNew >= fk) {
        this->success_ = false;
        return T(0.0);
    }
    this->success_ = true;
    return alpha;
}


template <typename T, typename F>
inline int LineSearch<T, F>::numFuncEvaluations() const {
    return funcNum_;
}

template <typename T, typename F>
inline bool LineSearch<T, F>::isSuccess() const {
    return success_;
}

#endif
