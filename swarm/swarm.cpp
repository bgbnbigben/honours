#include <swarm/swarm.h>
#include <swarm/particle.h>
#include <utilities/function.h>
#include <utilities/rrect.h>
#include <utilities/vector_ops.h>
#include <algorithm>

namespace {
    bool compare(Particle* a, Particle* b) { return a->getVal() < b->getVal();}
};

Swarm::Swarm(Function<double>* f, int n, int dim) : f_(f), number_(n), dimension_(dim), done_(false), same_(0), iterations_(0) {
    this->rtree_ = new RTree<double, double, DIMS>;
    this->particles_.resize(n);
    std::generate(this->particles_.begin(), this->particles_.end(), [=] () {
        return new Particle(Swarm::left_window, Swarm::right_window, dim);
    });
    this->setBests_();
}

void Swarm::setBests_() {
    double best = this->particles_[0]->getVal();
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [this] (Particle* p) {
            p->setVal((* this->f_)(p->pos()));
            RRect r(p->pos());
            this->rtree_->Insert(r.min, r.max, p->getVal());
    });
    std::sort(this->particles_.begin(), this->particles_.end(), ::compare);
    if (this->particles_[0]->getVal() > best)
        this->same_++;
    if (this->same_ > 100)
        this->done_ = true;
}

void Swarm::dance() {
    this->iterations_++;
    auto best = this->bestX();
    std::vector<double> centre(this->dimension_);
    # pragma omp parallel
    {
        std::for_each(this->particles_.begin(), this->particles_.end(),
            [&] (Particle*& i) {
                i->step(best, this->c1, this->c2, this->momentum);
                RRect r(i->pos());
                centre += i->pos();
                if (this->rtree_->Search(r.min, r.max, NULL, NULL))
                    i->reload(i->pos()+.2, i->vel()+.2);

        });
        centre /= double(this->number_);
    }
    this->setBests_();
    RRect r(centre);
    if (this->rtree_->Search(r.min, r.max, NULL, NULL) > this->number_ / 2) {
        this->done_ = true;
    }
}

double Swarm::bestVal() {
    return this->particles_[0]->getVal();
}

const std::vector<double>& Swarm::bestX() {
    return this->particles_[0]->pos();
}
