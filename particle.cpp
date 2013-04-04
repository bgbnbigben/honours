#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>
#include "particle.h"
#include "swarm.h"

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
    for (auto i = 0u; i < this->velocity_.size(); i++) {
        this->velocity_[i] = momentum * (this->velocity_[i] +
            c1 * ::getRandom() * (this->bestPosition_[i] - this->position_[i]) +
            c2 * ::getRandom() * (direction[i] - this->position_[i]));
        this->position_[i] += this->velocity_[i];
    }
}

int main() {
    Swarm swarm(Function(), 200, 50);
    double bestF = swarm.bestVal();
    int same = 0;
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(10);
    while (same < 10) {
        swarm.dance();
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
    std::cout << ")" << std::endl;
    return 0;
}
