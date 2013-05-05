#ifndef SWARM_H
#define SWARM_H
#include <vector>
#include <numeric>
#include "RTree.h"

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
        // TODO Refactor this properly. All the templates.
        int dimension_;
        const static int DIMS = 50;
        // <ID, data type, numDims>
        // Use the id as the function value maybe?
        RTree<double, double, DIMS>* rtree_;

    public:
        std::vector<Particle*> particles_;
    private:

        void setBests_();

    public:
        Swarm(Function, int, int);

        void dance();
        double bestVal();
        const std::vector<double>& bestX();
        const int numBlockedOff() { return rtree_->Count(); }

        static constexpr double left_window = -30.0;
        static constexpr double right_window = 30.0;
        static constexpr double c1 = 2.8;
        static constexpr double c2 = 1.3;
        static constexpr double momentum = 0.729844;

};

#endif
