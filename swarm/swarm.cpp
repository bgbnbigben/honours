#include <swarm/swarm.h>
#include <swarm/particle.h>
#include <utilities/function.h>
#include <utilities/vector_ops.h>
#include <utilities/RTreeUtilities.h>
#include <algorithm>

const double hyperrect_size = 0.5;

namespace {
    bool compare(Particle* a, Particle* b) { return a->getVal() < b->getVal();}
    bool boundCompare(const Bound<double> a, const Bound<double> b) { return a.variable < b.variable; }
};

Swarm::Swarm(Function<double>* f, int n, int dim, std::vector<Bound<double>> bounds) : f_(f), number_(n), dimension_(dim), done_(false), same_(0), iterations_(0), bounds_(bounds) {
    // The tree's ID doesn't matter, since we throw it away when it runs out of
    // scope. There's no need to reload it.
    SpatialIndex::id_type treeID;
    this->memoryStorage_ = SpatialIndex::StorageManager::createNewMemoryStorageManager();
    this->rtree_ = SpatialIndex::RTree::createNewRTree(*this->memoryStorage_, .7, 100, 100, dim, SpatialIndex::RTree::RV_RSTAR, treeID);
    this->particles_.resize(n);
    std::sort(this->bounds_.begin(), this->bounds_.end(), ::boundCompare);
    # pragma omp parallel
    {
        # pragma omp for
        for (int i = 0; i < this->particles_.size(); i++) {
            this->particles_[i] = new Particle(this->bounds_[i].lower, this->bounds_[i].upper, dim);
        }
    }
    this->setBests_();
}

Swarm::~Swarm() {
    delete this->rtree_;
    delete this->memoryStorage_;
}

void Swarm::setBests_() {
    double best = this->particles_[0]->getVal();
    std::vector<double> velocity(this->dimension_);
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [this, &velocity] (Particle* p) {
            p->setVal((* this->f_)(p->pos()));
            SpatialIndex::Region r(reinterpret_cast<const double*>(&(p->pos() - hyperrect_size)[0]), reinterpret_cast<const double*>(&(p->pos() + hyperrect_size)[0]), this->dimension_);
            velocity += p->vel();
            // We don't care to associate any data with this region
            this->rtree_->insertData(0, 0, r, 0);
    });
    std::sort(this->particles_.begin(), this->particles_.end(), ::compare);
    if (this->particles_[0]->getVal() > best)
        this->same_++;
    if (this->same_ > 100)
        this->done_ = true;
    if (norm(velocity) < 1.0)
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
            i->step(best, this->c1, this->c2, this->momentum, this->bounds_);
            centre += i->pos();
            IntersectingVisitor v;
            while (v.intersections()) {
                v.prepareForQuery();
                SpatialIndex::Region r(reinterpret_cast<const double*>(&i->pos()[0]), reinterpret_cast<const double*>(&i->pos()[0]), this->dimension_);
                this->rtree_->intersectsWithQuery(r, v);
                if (v.intersections()) {
                    i->reload(i->pos() + .2, i->vel() + .2);
                }
            }
        });
    }
    centre /= double(this->number_);
    this->setBests_();
    SpatialIndex::Region r(reinterpret_cast<const double*>(&(centre - hyperrect_size)[0]), reinterpret_cast<const double*>(&(centre + hyperrect_size)[0]), this->dimension_);
    IntersectingVisitor v;
    this->rtree_->intersectsWithQuery(r, v);
    if (v.intersections() > this->number_ / 2) {
        this->done_ = true;
    }
}

double Swarm::bestVal() {
    return this->particles_[0]->getVal();
}

const std::vector<double>& Swarm::bestX() {
    return this->particles_[0]->pos();
}
