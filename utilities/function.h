#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <algorithm>
#include <utilities/vector_ops.h>
#include <random>
#include <cassert>
#include <cutest.h>

template <class T>
class Function {
    protected:
        unsigned calls_;
        unsigned grads_;

    public:
        Function() : calls_(0), grads_(0) {};
        virtual T operator() (const std::vector<T>&) = 0;
        virtual std::vector<T> grad(const std::vector<T>&) = 0;
        inline unsigned numCalls() const { return this->calls_; }
        inline unsigned numGrads() const { return this->grads_; }
};

template <class T>
class Paraboloid : public Function<T> {
    public:
        Paraboloid() : Function<T>() {};

        T operator() (const std::vector<T> &x) {
            this->calls_++;
            return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        }

        std::vector<T> grad(const std::vector<T> &x) {
            this->grads_++;
            return T(2)*x;
        }
};

template <class T>
class Rosenbrock : public Function<T> {
    public:
        Rosenbrock() : Function<T>() {};
        T operator()(const std::vector<T> &x) {
            this->calls_++;
            T s = 0.0;
            for (int i = 0; i < x.size(); i++) {
                s += (1-x[i])*(1-x[i]);
                if (i < x.size() - 1)
                    s += 100 * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);
            }   
            return s;
        }   

        std::vector<T> grad(const std::vector<T> &x) {
            auto n = x.size();
            std::vector<T> df(n);
            std::vector<std::vector<T>> v(n, std::vector<T>(n, 0));
            for (int i = 0; i < n; i++) v[i][i] = .0000000000001;
            std::vector<T> delp(n);
            std::vector<T> delm(n);
            auto newx = x;
//            std::random_device rd;
//            std::mt19937 gen(rd());
//            std::uniform_real_distribution<> dis(-1, 1);
//            for (int i = 0; i < n; i++) {
//                newx[i] += dis(gen);
//            }
            auto xb = newx;
            auto fb = this->operator()(newx);
            auto fc = fb;
            int fdiff=0;
            for (int j = 0; j < n; j++) {
                auto xp = newx + v[j]; auto xm = newx - v[j]; auto fp = this->operator()(xp); delp[j] = fp - fc;
                if (fp < fb) fb = fp, xb = xp;
                if (fdiff == 0) { auto fm = this->operator()(xm); delm[j] = fc - fm;
                    if (fm < fb) fb = fm, xb = xm;
                }
            }
            auto invtran = inverse(transpose(v));
            if (fdiff == 1) {
                //Ax = b => x = A\b
                auto blah = invtran*delp;
                df = blah;
            } else
                df = .5 * (invtran*delp + invtran*delm);
            //this->grads_++;
            std::vector<T> df2(n);
            for (int i = 0; i < x.size(); i++) {
                if (i < x.size() - 1)
                    df2[i] += -400.0*x[i] * (x[i+1] - x[i] * x[i]) - 2.0*(1.0-x[i]);
                if (i > 0)
                    df2[i] += 200.0 * (x[i] - x[i-1]*x[i-1]);
            }
            return df2;
            std::cout << norm(df - df2) << std::endl;
            assert(norm(df - df2) <= 1);
            return df;
        }   
};

template <class T>
class CUTEst : public Function<T> {
    private:
        int nVar_;
        int nConstraints_;
        T *x_l, *x_u;

    protected:
        bool constrained_;
    public:
        CUTEst() : Function<T>(), constrained_(false) {
            char* fname = "tmp/OUTSDIF.d"; /* CUTEst data file */
            int funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
            int iout  = 6;         /* FORTRAN unit number for error output */
            int io_buffer = 11;    /* Internal input/output buffer */
            int ierr;              /* Exit flag for various calls */

            double *x; /* position, lower bound, upper bound */
            double *y = NULL, *y_l = NULL, *y_u = NULL;
            bool *equatn = NULL, *linear = NULL;
            int eqn_order = 0, l_order = 0, v_order = 0;

            /* Open problem description file OUTSDIF.d */
            ierr = 0;
            FORTRAN_open(&funit, fname, &ierr);
            if (ierr) {
                /* work it out later */
            }
            /* Determine problem size */
            CUTEST_cdimen(&ierr, &funit, &this->nVar_, &this->nConstraints_);
            if (ierr) {
                /* work it out later */
            }
            /* Determine whether to call constrained or unconstrained tools */
            if (this->nConstraints_) this->constrained_ = true;
            x = new double[this->nVar_];
            x_l = new T[this->nVar_];
            x_u = new T[this->nVar_];
            if (this->constrained_) {
                equatn = new bool[this->nConstraints_+1];
                linear = new bool[this->nConstraints_+1];
                y = new double[this->nConstraints_+1];
                y_l = new double[this->nConstraints_+1];
                y_u = new double[this->nConstraints_+1];
                CUTEST_csetup(&ierr, &funit, &iout, &io_buffer,
                        &this->nVar_, &this->nConstraints_,
                        x, x_l, x_u, y, y_l, y_u,
                        equatn, linear, &eqn_order, &l_order, &v_order);
                if (ierr) {
                    /* work it out later */
                }
            } else {
                equatn = new bool[1];
                linear = new bool[1];
                y_l = new double[1];
                y_u = new double[1];
                CUTEST_usetup(&ierr, &funit, &iout, &io_buffer,
                        &this->nVar_, x, x_l, x_u);
                if (ierr) {
                    /* work it out later */
                }
            }
            FORTRAN_close(&funit, &ierr);
            delete[] x;
            delete[] y;
            delete[] y_l;
            delete[] y_u;
            delete[] equatn;
            delete[] linear;
        }

        T operator() (const std::vector<T> &x) {
            this->calls_++;
            T* grad = NULL;
            T obj;
            int ierr;
            if (this->constrained_) {
                logical a = false;
                CUTEST_cofg(&ierr, &this->nVar_, const_cast<T*>(&x.front()), &obj, grad, &a);
            } else {
                logical a = false;
                CUTEST_uofg(&ierr, &this->nVar_, const_cast<T*>(&x.front()), &obj, grad, &a);
            }
            return obj;
        }

        std::vector<T> grad(const std::vector<T> &x) {
            this->grads_++;
            std::vector<T> grad(this->nVar_);
            T obj;
            int ierr;
            if (this->constrained_) {
                logical a = false;
                CUTEST_cofg(&ierr, &this->nVar_, const_cast<T*>(&x.front()), &obj, &grad.front(), &a);
            } else {
                logical a = false;
                CUTEST_uofg(&ierr, &this->nVar_, const_cast<T*>(&x.front()), &obj, &grad.front(), &a);
            }
            return grad;
        }

        ~CUTEst() {
            delete[] x_l;
            delete[] x_u;
        }

        const int getDim() { return this->nVar_; }

        const std::vector<T> getLowerBounds() {
            return std::vector<T>(x_l, x_l + this->nVar_);
        }

        const std::vector<T> getUpperBounds() {
            return std::vector<T>(x_u, x_u + this->nVar_);
        }
};
#endif
