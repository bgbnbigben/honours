#ifndef SWARM_H
#define SWARM_H
#include <vector>
#include <numeric>

class Particle;

class Function {
    public:
        double operator() (std::vector<double> x) {
            return std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        }
};

class Swarm {
        Function f_;
        int number_;
        int dimension_;
        std::vector<Particle*> particles_;

        void setBests_();

    public:
        Swarm(Function, int, int);

        void dance();
        double bestVal();
        const std::vector<double>& bestX();

        static const double left_window = -30.0;
        static const double right_window = 30.0;
        static const double c1 = 2.8;
        static const double c2 = 1.3;
        static const double momentum = 0.729844;
};

#endif
