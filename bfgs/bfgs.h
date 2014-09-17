#ifndef BFGS_H
#define BFGS_H

#include <utilities/bound.h>
#include <utilities/function.h>
#include <bfgs/lbfgsb_prototypes.h>
#include <cstring>

template <typename T>
class BFGS { 
    public:

    void optimize(Function<T>&, std::vector<T>&, std::vector<Bound<T>> b, T decreaseFactor=1e1,
                  T projectionTol=1e-5);

        std::vector<T> getOptValue() const;
        T getFuncMin() const;
        int numIterations() const;
        inline bool isSuccess() const { return success_; }

    private:

        T fMin_;
        std::vector<T> xOpt_;
        std::vector<T> gradNorm_;
        int iters_;
        bool success_;
};

/**
 */
template <typename T>
void BFGS<T>::optimize(Function<T> &func, std::vector<T> &x0, std::vector<Bound<T>> bounds, T decreaseFactor, T projectionTol) {
    //char task[60 + 1 + 1], csave[60 + 1 + 1];
    char task[60];
    char csave[60];
    for (int i = 0; i < 60; i++) csave[i] = ' ';
    strncpy(task, "START", 60);
    for (int i = 5; i < 60; i++) task[i] = ' ';

    char lsave[4] = {0, 0, 0, 0};
    int n = x0.size(), m = 17, iprint = 1; 
    int isave[44];

    double *x = new double[n];
    std::copy(x0.begin(), x0.end(), x);

    int *nbd = new int[n];
    for (int i = 0; i < n; i++) nbd[i] = 0;
    double *l = new double[n]; 
    double *u = new double[n];
    for (auto bound: bounds) {
        nbd[bound.variable] = bound.type;
        if (bound.type == Bound<T>::types::LOWER || bound.type == Bound<T>::types::BOTH)
            l[bound.variable] = bound.lower;
        if (bound.type == Bound<T>::types::UPPER || bound.type == Bound<T>::types::BOTH)
            u[bound.variable] = bound.upper;
    }

    int *iwa = new int[3*n];

    double f, dsave[29];
    double *g = new double[n];
    double *wa = new double[2*m*n+5*n+11*m*m+8*m];

    while (true) {
        setulb_(&n, &m, x, l, u, nbd, &f, g, &decreaseFactor, &projectionTol, wa, iwa, task, &iprint, csave, lsave, isave, dsave, 60, 60, 4);
        std::copy(x, x+n, x0.begin());

        if (task[0] == 'F' && task[1] == 'G') {
            /* Wants f and g for x */
            f = func(x0);
            auto temp = func.grad(x0);
            std::copy(temp.begin(), temp.end(), g);
        } else if (task[0] == 'N' && 
                   task[1] == 'E' &&
                   task[2] == 'W' &&
                   task[3] == '_' &&
                   task[4] == 'X') {
            /* Any stopping criteria can go here; also note that you must set
             * task[] to "STOP"
             */
        } else {
            xOpt_.resize(n);
            std::copy(x, x+n, xOpt_.begin());
            fMin_ = f;
            iters_ = isave[30];
            if (task[0] == 'C' && 
                task[1] == 'O' && 
                task[2] == 'N' && 
                task[3] == 'V')
                success_ = true;
            else success_ = false;
            task[59] = 0;
            //std::cout << task << std::endl;

            delete[] x;
            delete[] nbd;
            delete[] l;
            delete[] u;
            delete[] iwa;
            delete[] g;
            delete[] wa;
            return;
        }
    }
}

/**
 * Get the optimum point.
 */
template <typename T>
inline std::vector<T> BFGS<T>::getOptValue() const {
    return xOpt_;
}


/**
 * Get the minimum value of objective function.
 */
template <typename T>
inline T BFGS<T>::getFuncMin() const {
    return fMin_;
}


/**
 * Get the iteration number.
 */
template <typename T>
inline int BFGS<T>::numIterations() const {
    return iters_;
}

#endif
// BFGS_H
