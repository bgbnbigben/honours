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
Particle<T>::Particle(const std::vector<VariableContainer>& position, int dim) : 
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

template <typename T>
Variable Particle<T>::operator[] (int i) {
    return this->position_[i];
}

/* TODO: make it clamp to window or re-init or something */
template <typename T>
void Particle<T>::step(const std::vector<Variable>& direction, double c1, double c2, double momentum) {
    # pragma omp single
    {
        this->velocity_ = momentum * (this->velocity_ +
                c1 * ::getRandom() * (this->bestPosition_ - this->position_) +
                c2 * ::getRandom() * (direction - this->position_));
        this->position_ += this->velocity_;
    }
    this->clamp();
}

template <typename T>
void Particle<T>::clamp() {
    # pragma omp parallel
    {
        std::for_each(this->position_.begin(), this->position_.end(),
            [&] (VariableContainer& pos) {
                    if (pos.type() == VariableType::Continuous) {
                        pos.val(std::max((double)pos.rangeStart, (double)pos));
                        pos.val(std::min((double)pos.rangeEnd, (double)pos));
                    } else {
                        pos.val(std::max((long long)pos.rangeStart, (long long)pos));
                        pos.val(std::min((long long)pos.rangeEnd, (long long)pos));
                        // Also force it to be at an integer point
                        pos.val((long long)((double)((long long)pos) + .5));
                    }
            });
    }
}

template <typename T>
void Particle<T>::reload(const std::vector<Variable>& p, const std::vector<Variable>& v) {
    this->velocity_ = v;
    this->position_ = p;
    this->clamp();
}

template <typename T>
VariableType Particle<T>::getRandomInWindow_(int i) {
    if (this->position_[i].getType == VariableType::Continuous)
        return ::getRandom() * (this->position_[i].rangeEnd.d - this->position_[i].rangeStart.d)
            + this->position_[i].rangeStart.d;
        return ::getRandom() * (this->position_[i].rangeEnd.ll - this->position_[i].rangeStart.ll)
            + this->position_[i].rangeStart.ll;
}
