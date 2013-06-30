#include <swarm/swarm.h>
#include <swarm/particle.h>
#include <utilities/function.h>
#include <utilities/rrect.h>
#include <algorithm>

namespace {
    bool compare(Particle* a, Particle* b) { return a->getVal() < b->getVal();}
};

Swarm::Swarm(Function<double>* f, int n, int dim) : f_(f), number_(n), dimension_(dim) {
    this->rtree_ = new RTree<double, double, DIMS>;
    this->particles_.resize(n);
    std::generate(this->particles_.begin(), this->particles_.end(), [=] () {
        return new Particle(Swarm::left_window, Swarm::right_window, dim);
    });
    this->setBests_();
}

void Swarm::setBests_() {
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [this] (Particle* p) {
            p->setVal((* this->f_)(p->pos()));
            RRect r(p->pos());
            this->rtree_->Insert(r.min, r.max, p->getVal());
    });
    std::sort(this->particles_.begin(), this->particles_.end(), ::compare);
}

#include <iostream>

void Swarm::dance() {
    std::vector<double> best = this->bestX();
    # pragma omp parallel
    {
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [&] (Particle*& i) {
            i->step(best, this->c1, this->c2, this->momentum);
            RRect r(i->pos());
            if (this->rtree_->Search(r.min, r.max, NULL, NULL)) {
                std::cout << "Duplicate! This rect intersects " << this->rtree_->Search(r.min, r.max, NULL, NULL) << " rectangles" << std::endl;
                delete i;
                i = new Particle(Swarm::left_window, Swarm::right_window, this->dimension_);
            }
    });
    }
    this->setBests_();
}

double Swarm::bestVal() {
    return this->particles_[0]->getVal();
}

const std::vector<double>& Swarm::bestX() {
    return this->particles_[0]->pos();
}
