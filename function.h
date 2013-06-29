#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <algorithm>

template <class type>
class Function {
    public:
        virtual type operator() (const std::vector<type>&) = 0;
        virtual std::vector<type> grad(const std::vector<type>&) = 0;
};

template <class type>
class Paraboloid : public Function<type> {
    public:
        type operator() (const std::vector<type> &x) {
            return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        }

        std::vector<type> grad(const std::vector<type> &x) {
            std::vector<type> r(x.size());
            std::transform(x.begin(), x.end(), r, [](type& i) {
                return type(2)*i;
            });
            return r;
        }
};

template <class type>
class Rosenbrock : public Function<type> {
    public:
        type operator()(const std::vector<type> &x) {
            type s = 0.0;
            for (int i = 0; i < x.size(); i++) {
                s += (1-x[i])*(1-x[i]);
                if (i < x.size() - 1)
                    s += 100 * (x[i+1] - x[i]*x[i]) * (x[i+1] - x[i]*x[i]);
            }   
            return s;
        }   

        std::vector<type> grad(const std::vector<type> &x) {
            std::vector<type> df(x.size());
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
