#ifndef SWARM_H
#define SWARM_H
#include <vector>
#include <numeric>
#include <utilities/RTree.h>
#include <utilities/function.h>

class Particle;

class Swarm {
        Function<double>* f_;
        unsigned number_;
        // TODO Refactor this properly. All the templates.
        int dimension_;
        const static int DIMS = 10;
        // <ID, data type, numDims>
        // Use the id as the function value maybe?
        RTree<double, double, DIMS>* rtree_;
        bool done_;
        int same_;
        unsigned iterations_;

    public:
        std::vector<Particle*> particles_;
    private:

        void setBests_();

    public:
        Swarm(Function<double>*, int, int);

        void dance();
        double bestVal();
        const std::vector<double>& bestX();
        bool done() const { return this->done_; }

        int numIterations() const { return this->iterations_; }
        int numParticles() const { return this->number_; }

        static constexpr double left_window = -30.0;
        static constexpr double right_window = 30.0;
        static constexpr double c1 = 2.8;
        static constexpr double c2 = 1.3;
        static constexpr double momentum = 0.729844;

};

#endif
