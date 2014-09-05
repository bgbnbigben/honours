#include <algorithm>
#include <omp.h>
#include <random>
#include <cmath>
#include <swarm/particle.h>
#include <swarm/neldermead.h>
#include <swarm/swarm.h>
#include <utilities/vector_ops.h>
#include <utilities/bound.h>

namespace {
    std::random_device rand;
    std::uniform_real_distribution<double> uniform;

    double getRandom() {
        return uniform(rand);
    }
}

Particle::Particle(const std::vector<Bound<double>>& bounds, int dim) : 
        bounds_(bounds) {
    assert(bounds.size() == dim);
    this->position_.resize(dim);
    this->velocity_.resize(dim);
    this->bestPosition_.resize(dim);
    # pragma omp parallel for
        for (int i = 0; i < dim; i++)
            this->position_[i] = this->getRandomInWindow_(i);
    # pragma omp parallel for
        for (int i = 0; i < dim; i++)
            this->velocity_[i] = this->getRandomInWindow_(i);
    this->bestPosition_ = this->position_;
    this->cost_ = 1e100;
    this->clamp();
}

double& Particle::operator[] (int i) {
    return this->position_[i];
}

/* TODO: make it clamp to window or re-init or something */
void Particle::step(const std::vector<double>& direction, double c1, double c2, double momentum) {
    # pragma omp single
    {
        this->velocity_ = momentum * (this->velocity_ +
                c1 * ::getRandom() * (this->bestPosition_ - this->position_) +
                c2 * ::getRandom() * (direction - this->position_));
        this->position_ += this->velocity_;
    }
    this->clamp();
}

void Particle::clamp() {
    # pragma omp parallel
    {
        std::for_each(this->bounds_.begin(), this->bounds_.end(),
            [&] (Bound<double> bound) {
                if (bound.type == Bound<double>::types::LOWER ||
                    bound.type == Bound<double>::types::BOTH) 
                    this->position_[bound.variable] = std::max(this->position_[bound.variable], bound.lower);
                if (bound.type == Bound<double>::types::UPPER ||
                    bound.type == Bound<double>::types::BOTH)
                    this->position_[bound.variable] = std::min(this->position_[bound.variable], bound.upper);
            });
    }
}

void Particle::reload(const std::vector<double>& p, const std::vector<double>& v) {
    this->velocity_ = v;
    this->position_ = p;
    this->clamp();
}

double Particle::getRandomInWindow_(int i) {
    return ::getRandom() * (this->bounds_[i].upper - this->bounds_[i].lower)
        + this->bounds_[i].lower;
}
