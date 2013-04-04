#include "swarm.h"
#include "particle.h"
#include <algorithm>

Swarm::Swarm(Function f, int n, int dim) : f_(f), number_(n), dimension_(dim) {
    this->particles_.resize(n);
    std::generate(this->particles_.begin(), this->particles_.end(), [=] () {
        return new Particle(Swarm::left_window, Swarm::right_window, dim);
    });
    this->setBests_();
}

bool compare(Particle* a, Particle* b) { return a->getVal() < b->getVal();}

void Swarm::setBests_() {
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [this] (Particle* p) {
            p->setVal(this->f_(p->pos()));
    });
    std::sort(this->particles_.begin(), this->particles_.end(), compare);
}

void Swarm::dance() {
    std::vector<double> best = this->bestX();
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [&] (Particle* i) {
            i->step(best, this->c1, this->c2, this->momentum);
    });
    this->setBests_();
}

double Swarm::bestVal() {
    return this->particles_[0]->getVal();
}

const std::vector<double>& Swarm::bestX() {
    return this->particles_[0]->pos();
}
