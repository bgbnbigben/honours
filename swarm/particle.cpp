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

template <typename T>
Particle::Particle(const std::vector<VariableContainer>& position, int dim) : 
        position_(position) {
    assert(position.size() == dim);
    this->velocity_.resize(dim);
    this->bestPosition_.resize(dim);
    //# pragma omp parallel for
    //    for (int i = 0; i < dim; i++)
    //        this->position_[i] = this->getRandomInWindow_(i);
    # pragma omp parallel for
        for (int i = 0; i < dim; i++) {
            auto r = this->getRandomInWindow_(i);
            this->velocity_[i] = this->position_[i].getType == VariableType::Continuous ? r.d : r.ll;
            this->velocity_[i].setType(this->position_[i].getType());
        }
    this->bestPosition_ = this->position_;
    this->cost_ = (T)1e100;
    this->clamp();
}

Variable& Particle::operator[] (int i) {
    return this->position_[i];
}

/* TODO: make it clamp to window or re-init or something */
void Particle::step(const std::vector<Variable>& direction, double c1, double c2, double momentum) {
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
        std::for_each(this->position_.begin(), this->position_.end(),
            [&] (VariableContainer& pos) {
                    // TODO make sure it's not fucking up and calling the wrong
                    // cast operator. It probably always uses double.
                    pos.val(std::max(pos.rangeStart, pos));
                    pos.val(std::min(pos.rangeEnd, pos));
            });
    }
}

void Particle::reload(const std::vector<Variable>& p, const std::vector<Variable>& v) {
    this->velocity_ = v;
    this->position_ = p;
    this->clamp();
}

VariableType Particle::getRandomInWindow_(int i) {
    if (this->position_[i].getType == VariableType::Continuous)
        return ::getRandom() * (this->position_[i].rangeEnd.d - this->position_[i].rangeStart.d)
            + this->position_[i].rangeStart.d;
        return ::getRandom() * (this->position_[i].rangeEnd.ll - this->position_[i].rangeStart.ll)
            + this->position_[i].rangeStart.ll;
}
