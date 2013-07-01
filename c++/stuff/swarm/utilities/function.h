#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <algorithm>

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
            this->grads_++;
            std::vector<T> df(x.size());
            for (int i = 0; i < x.size(); i++) {
                if (i < x.size() - 1)
                    df[i] = -400*x[i] * (x[i+1] - x[i] * x[i]) -2 * (1-x[i]);
                if (i > 0)
                    df[i] += 200 * (x[i] - x[i-1]*x[i-1]);
            }   
            return df; 
        }   
};
#endif
