#include <swarm/neldermead.h>
#include <swarm/swarm.h>
#include <utilities/vector_ops.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace {
    class NumIterationsExceeded : public std::length_error {
        public:
            NumIterationsExceeded(const char* w) : std::length_error(w) {};
    };
};


NelderMead::NelderMead(int dim, double tol, int iter) : dimension_(dim), termination_tolerance_(tol), max_num_iterations_(iter), iteration_count_(0) {
    this->points_.reserve(this->dimension_ + 1);
}

bool NelderMead::operator()(const std::vector<double>& a, const std::vector<double>& b) {
    return lookup(a) < lookup(b);
}

bool NelderMead::done_() {
    if (this->max_num_iterations_ && this->iteration_count_ == this->max_num_iterations_) throw NumIterationsExceeded("Number of iterations has been exceeded!");
    if (this->points_.size() < this->dimension_) return false;
    for (unsigned i = 0; i < this->dimension_ + 1; i++)
        for (unsigned j = i + 1; j < this->dimension_ + 1; j++)
            if (std::inner_product(this->points_[i].begin(), this->points_[i].end(), this->points_[j].begin(), 0.0, std::plus<double>(), [] (double a, double b) { return (a-b)*(a-b);}) > this->termination_tolerance_ * this->termination_tolerance_) return false;

    return true;
}

double NelderMead::lookup(const std::vector<double>& point) const {
    auto pos = this->cache_.find(point);
    if (pos != this->cache_.end())
        return pos->second;
    throw point;
}

void NelderMead::addSimplexPoint(std::vector<double> point) {
    /* Silently fail */
    if (this->points_.size() < this->dimension_ + 1)
        this->points_.push_back(point);
}

std::vector<double> NelderMead::calc_centre_(bool end) {
    std::vector<double> centre(this->dimension_, 0);
    for (unsigned i = 0; i < this->dimension_; i++) {
        for (unsigned j = 0; j < this->points_.size(); j++) {
            // Sum stupidly so we can do coordinate-wise
            centre[i] += this->points_[j][i];
        }
    }
    std::for_each(centre.begin(), centre.end(), [&](double& pt) {
        pt /= (double)(this->dimension_ - (end ? 0.0 : 1.0));
    });
    return centre;
}

std::vector<double> NelderMead::nextIteration(std::vector<double> point, double f) {
    this->cache_[point] = f;
    try { 
        if (this->points_.size() >= this->dimension_ + 1) {
            while (!this->done_()) {
                this->iteration_count_++;
                std::sort(this->points_.begin(), this->points_.end(), *this);
                std::vector<double> centre = calc_centre_(false);

                std::vector<double> best(this->points_[0]);
                std::vector<double> worst(this->points_[this->dimension_]);
                std::vector<double> second_worst(this->points_[this->dimension_ - 1]);

                std::vector<double> reflected = centre + NelderMead::alpha * (centre - worst);
                if (lookup(reflected) < lookup(second_worst) && lookup(reflected) > lookup(best)) {
                    this->points_[this->dimension_] = reflected;
                } else if (lookup(reflected) < lookup(best)) {
                    std::vector<double> expanded = centre + NelderMead::gamma * (centre - worst);
                    this->points_[this->dimension_] = (lookup(expanded) > lookup(reflected)) ? expanded : reflected;
                } else {
                    std::vector<double> contracted = centre - NelderMead::rho * (centre - worst);
                    if (lookup(contracted) < lookup(worst)) {
                        this->points_[this->dimension_] = contracted;
                    } else {
                        std::for_each(this->points_.begin(), this->points_.end(), [&](std::vector<double>& point) {
                            point = best + NelderMead::sigma * (point - best);
                        });
                    }
                }
            }

            return calc_centre_(true);
        } else {
            std::cerr << "Error condition" << std::endl;
            return std::vector<double>(this->dimension_, -1000000.0);
        }
    } catch (std::vector<double> v) {
        // lookup throws a vector if it doesn't exist in the cache
        return v;
    }
}

std::vector<double> NelderMead::drive(Function<double>* f) {
    std::vector<double> curr = this->points_[0];
    try {
        while (!this->done_()) {
            curr = this->nextIteration(curr, (*f)(curr));
            //print(curr)
        }  
    } catch (NumIterationsExceeded& n) {
        curr = this->calc_centre_(true);
    }
    return curr;
}
