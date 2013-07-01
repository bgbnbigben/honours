#ifndef BFGS_H
#define BFGS_H

#include <utilities/matrix.h>
#include <bfgs/linesearch.h>
#include <limits>

template <typename T, typename F>
class BFGS : public LineSearch<T, F> {
    public:

        BFGS(): LineSearch<T, F>() {};
        ~BFGS() {};

        void optimize(F &func, std::vector<T> &x0, T tol=T(1.0e-6), int maxItr=100);

        std::vector<T> getOptValue() const;
        std::vector<T> getGradNorm() const;
        T getFuncMin() const;
        int numIterations() const;

    private:

        T fMin_;
        std::vector<T> xOpt_;
        std::vector<T> gradNorm_;

        const T EPS_ = std::numeric_limits<T>::epsilon();

};

/**
 * Finding the optimal solution. The default tolerance error and maximum
 * iteratin number are "tol=1.0e-6" and "maxItr=100", respectively.
 */
template <typename T, typename F>
void BFGS<T, F>::optimize(F &func, std::vector<T> &x0, T tol, int maxItr) {
    // initialize parameters.
    unsigned k = 0, cnt = 0, N = x0.size();

    T ys, yHy, alpha;
    std::vector<T> d(N),
                   s(N),
                   y(N),
                   v(N),
                   Hy(N),
                   gPrev(N);
    Matrix<T> H = eye<T>(N);

    std::vector<T> x(x0);
    T fx = func(x);
    this->funcNum_++;
    std::vector<T> gnorm(maxItr);
    std::vector<T> g = func.grad(x);
    gnorm[k++]= norm(g);

    while ((gnorm[k-1] > tol) && (k < maxItr)) {
        // descent direction
        d = - H * g;

        // one dimension searching
        alpha = this->getStep(func, x, d);

        // check flag for restart
        if (!this->success_)
            // Test if the norm of (H - I)
            if (norm(H - eye<T>(N)) < this->EPS_)
                break;
            else {
                H = eye<T>(N);
                cnt++;
                if (cnt == maxItr)
                    break;
            }
        else {
            // update
            s = alpha * d;
            x += s;
            fx = func(x);
            this->funcNum_++;
            gPrev = g;
            g = func.grad(x);
            y = g - gPrev;

            Hy = H * y;
            ys = dotProd(y, s);
            yHy = dotProd(y, Hy);

            // TODO: less bullshit / more fabs.
            if ((ys < this->EPS_) || (yHy < this->EPS_))
                H = eye<T>(N);
            else {
                v = sqrt(yHy) * (s/ys - Hy/yHy);
                H = H + multTr(s,s)/ys - multTr(Hy,Hy)/yHy + multTr(v,v);
            }
            gnorm[k++] = norm(g);
        }
    }

    xOpt_ = x;
    fMin_ = fx;
    gradNorm_.resize(k);
    for (unsigned i = 0; i < k; ++i)
        gradNorm_[i] = gnorm[i];

    if (gradNorm_[k-1] > tol)
        this->success_ = false;
}


/**
 * Get the optimum point.
 */
template <typename T, typename F>
inline std::vector<T> BFGS<T, F>::getOptValue() const {
    return xOpt_;
}


/**
 * Get the norm of gradient in each iteration.
 */
template <typename T, typename F>
inline std::vector<T> BFGS<T, F>::getGradNorm() const {
    return gradNorm_;
}


/**
 * Get the minimum value of objective function.
 */
template <typename T, typename F>
inline T BFGS<T, F>::getFuncMin() const {
    return fMin_;
}


/**
 * Get the iteration number.
 */
template <typename T, typename F>
inline int BFGS<T, F>::numIterations() const {
    return gradNorm_.size() - 1;
}

#endif
// BFGS_H
