#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "particle.h"
#include "swarm.h"
#include "neldermead.h"

namespace {
    std::random_device rand;
    std::uniform_real_distribution<double> uniform;

    double getRandom() {
        return uniform(rand);
    }

    double getRandomInWindow() {
        return ::getRandom() * (Swarm::right_window - Swarm::left_window)
            + Swarm::left_window;
    }

};

Particle::Particle(double left, double right, int dim) : 
        leftWindow_(left), rightWindow_(right) {
    this->position_.resize(dim);
    this->velocity_.resize(dim);
    this->bestPosition_.resize(dim);
    std::generate(this->position_.begin(),
            this->position_.end(),
            getRandomInWindow);
    std::generate(this->velocity_.begin(),
            this->velocity_.end(),
            getRandomInWindow);
    this->bestPosition_ = this->position_;
    this->cost_ = 1e100;
}

double& Particle::operator[] (int i) {
    return this->position_[i];
}

/* TODO: make it clamp to window or re-init or something */
void Particle::step(const std::vector<double>& direction, double c1, double c2, double momentum) {
    //if (getRandom() < Particle::mutation_prob) {
    //    NelderMead n(this->position_.size(), 0, 80);
    //    n.addSimplexPoint(this->position_);
    //    for (auto i = 0u; i < this->position_.size(); i++) {
    //        std::vector<double> perturbed = this->position_;
    //        perturbed[i] += (std::fabs(perturbed[i]) > 1e-9 ? 0.05 : 0.00025);
    //        n.addSimplexPoint(perturbed);
    //    }
    //    this->position_ = n.drive(new Function());
    //} else {
        for (auto i = 0u; i < this->velocity_.size(); i++) {
            this->velocity_[i] = momentum * (this->velocity_[i] +
                c1 * ::getRandom() * (this->bestPosition_[i] - this->position_[i]) +
                c2 * ::getRandom() * (direction[i] - this->position_[i]));
            this->position_[i] += this->velocity_[i];
        }
    //}
}

constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;

int main() {
    //NelderMead n(8, 0, 320);
    //n.addSimplexPoint({10, 10, 10, 10, 10, 10, 10, 10});
    //n.addSimplexPoint({-10, -10, -10, -10, -10, -10, -10, -10});
    //n.addSimplexPoint({20, 20, 20, 20, 20, 20, 20, 20});
    //n.addSimplexPoint({-20, -20, -20, -20, -20, -20, -20, -20});
    //n.addSimplexPoint({30, 30, 30, 30, 30, 30, 30, 30});
    //n.addSimplexPoint({-30, -30, -30, -30, -30, -30, -30, -30});
    //n.addSimplexPoint({40, 40, 40, 40, 40, 40, 40, 40});
    //n.addSimplexPoint({-40, -40, -40, -40, -40, -40, -40, -40});
    //n.addSimplexPoint({50, 50, 50, 50, 50, 50, 50, 50});
    //std::vector<double> soln = n.drive(new Function());
    //std::cout << (*(new Function()))(soln) << ":\t";
    //std::for_each(soln.begin(), soln.end(), [](double i) {
    //    std::cout << i << "\t";
    //});
    //std::cout << std::endl;

    Swarm swarm(Function(), 300, 50);
    double bestF = swarm.bestVal();
    int same = 0;
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(10);
    while (same < 10) {
        swarm.dance();
        std::cout << "Best in swarm: " << swarm.bestVal() << std::endl;
        std::cout << "Best so far: " << bestF << std::endl;
        if (bestF > swarm.bestVal()) {
            bestF = swarm.bestVal();
            same = 0;
        } else {
            same++;
        }
    }
    std::cout << bestF << ":\t( ";
    std::for_each(swarm.bestX().begin(), swarm.bestX().end(), [](double i) {
        std::cout << i << ",\t";
    });
    std::cout << ")\n" << std::endl;

    std::cout << "There are " << swarm.numBlockedOff() << " items in the tabu list" << std::endl;
    return 0;
}
