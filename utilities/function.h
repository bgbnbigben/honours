#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <algorithm>
#include <utilities/vector_ops.h>
#include <random>

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
            std::cout << norm(df - df2) << std::endl;
            assert(norm(df - df2) <= 1);
            return df;
        }   
};
#endif
