#include <algorithm>
#include <omp.h>
#include <random>
#include <cmath>
#include <swarm/particle.h>
#include <swarm/neldermead.h>
#include <swarm/swarm.h>
#include <utilities/vector_ops.h>


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
    //    NelderMead n(this->position_.size(), 0, 20);
    //    n.addSimplexPoint(this->position_);
    //    for (auto i = 0u; i < this->position_.size(); i++) {
    //        std::vector<double> perturbed = this->position_;
    //        perturbed[i] += (std::fabs(perturbed[i]) > 1e-9 ? 0.05 : 0.00025);
    //        n.addSimplexPoint(perturbed);
    //    }
    //    this->position_ = n.drive(new Function());
    //} else {
        # pragma omp single
        {
            this->velocity_ = momentum * (this->velocity_ +
                    c1 * ::getRandom() * (this->bestPosition_ - this->position_) +
                    c2 * ::getRandom() * (direction - this->position_));
            this->position_ += this->velocity_;
        }
    //}
}
