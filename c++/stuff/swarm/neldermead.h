#ifndef NELDERMEAD_H
#define NELDERMEAD_H
#include <vector>
#include <map>

class Function;

class NelderMead {
    private:
        unsigned dimension_;
        double termination_tolerance_;
        unsigned max_num_iterations_;
        unsigned iteration_count_;
        std::vector<std::vector<double>> points_;
        std::map<std::vector<double>, double> cache_;

        bool done_();
        double lookup(const std::vector<double>&) const;
        std::vector<double> calc_centre_(bool);

    public:
        NelderMead(int, double, int iter=10);
        bool operator()(const std::vector<double>&, const std::vector<double>&);
        void addSimplexPoint(std::vector<double>);
        std::vector<double> nextIteration(std::vector<double>, double); 
        std::vector<double> drive(Function*);

        static constexpr double alpha = 1;
        static constexpr double gamma = (1 + 2.0/50.0);
        static constexpr double rho   = (0.75 - 1.0/(2*50.0));
        static constexpr double sigma = -(1.0 - 1.0/50.0);
};

#endif
