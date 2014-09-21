#include <swarm/swarm.h>
#include <swarm/particle.h>
#include <utilities/function.h>
#include <utilities/vector_ops.h>
#include <utilities/RTreeUtilities.h>
#include <algorithm>
#include <unistd.h>
#include <iostream>
#include <fstream>

const double hyperrect_size = 0.5;

namespace {
    template <typename T>
    bool compare(Particle<T>* a, Particle<T>* b) { return a->getVal() < b->getVal();}

    template <typename T>
    bool inBounds(Particle<T>* p) {
        bool inBounds = true;
        for (unsigned i = 0; i < p->bounds_.size(); i++)
            inBounds = inBounds && p->bounds_[i].lower <= p->position_[i] && p->position_[i] <= p->bounds_[i].upper;
        return inBounds;
    }
};

template <typename T>
Swarm<T>::Swarm(Function<T>* f, int n, int dim, std::vector<VariableContainer> bounds) : f_(f), number_(n), dimension_(dim), done_(false), same_(0), iterations_(0) {
    // I probably will regret this.
    auto name = std::string("temp" + std::to_string(rand()));
    auto pageSize = sysconf(_SC_PAGESIZE);
    //std::cout << "The page size is " << pageSize << std::endl;
    // The tree's ID doesn't matter, since we throw it away when it runs out of
    // scope. There's no need to reload it.
    SpatialIndex::id_type treeID;
    //this->memoryStorage_ = SpatialIndex::StorageManager::createNewMemoryStorageManager();
    this->memoryStorage_ = SpatialIndex::StorageManager::createNewDiskStorageManager(name, pageSize);
    this->memoryBuffer_ = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*this->memoryStorage_, 10, false);
    this->rtree_ = SpatialIndex::RTree::createNewRTree(*this->memoryStorage_, .7, 100, 100, dim + (dim == 1), SpatialIndex::RTree::RV_RSTAR, treeID);
    this->particles_.resize(n);
    # pragma omp parallel
    {
        # pragma omp for
        for (int i = 0; i < this->particles_.size(); i++) {
            this->particles_[i] = new Particle<T>(bounds, dim);
        }
    }
    this->setBests_();
}

template <typename T>
Swarm<T>::~Swarm() {
    delete this->rtree_;
    delete this->memoryBuffer_;
    delete this->memoryStorage_;
}

template <typename T>
void Swarm<T>::setBests_() {
    T best = this->particles_[0]->getVal();
    std::vector<Variable> velocity(this->dimension_);
    std::for_each(this->particles_.begin(), this->particles_.end(),
        [this, &velocity] (Particle<T>* p) {
            p->setVal((* this->f_)(p->pos()));
            /* TODO use the DTRegion
            SpatialIndex::Region r;
            if (this->dimension_ == 1) {
                double lower[] = {(p->pos() - hyperrect_size)[0], 0};
                double upper[] = {(p->pos() + hyperrect_size)[0], 0};
                r = SpatialIndex::Region(lower, upper, 2);
            }
            else
                r = SpatialIndex::Region(reinterpret_cast<const double*>(&(p->pos() - hyperrect_size)[0]), reinterpret_cast<const double*>(&(p->pos() + hyperrect_size)[0]), this->dimension_);
            // We don't care to associate any data with this region
            this->rtree_->insertData(0, 0, r, 0);
            */
            velocity += p->vel();
    });
    std::sort(this->particles_.begin(), this->particles_.end(), ::compare);
    if (this->particles_[0]->getVal() > best)
        this->same_++;
    if (this->same_ > 100)
        this->done_ = true;
    // TODO use something meaningful for type VariableContainer
    //if (norm(velocity) < 1.0)
    //    this->done_ = true;
    if (this->iterations_ >= 1000)
        this->done_ = true;

}

template <typename T>
void Swarm<T>::dance() {
    this->iterations_++;
    auto best = this->bestX();
    std::vector<double> centre(this->dimension_);
    # pragma omp parallel
    {
        std::for_each(this->particles_.begin(), this->particles_.end(),
          [&] (Particle<T>*& i) {
            i->step(best, this->c1, this->c2, this->momentum);
            centre += i->pos();
            IntersectingVisitor v;
            while (v.intersections()) {
                v.prepareForQuery();
                SpatialIndex::Region r;
                if (this->dimension_ == 1) {
                    double lower[] = {i->pos()[0], 0};
                    double upper[] = {i->pos()[0], 0};
                    r = SpatialIndex::Region(lower, upper, 2);
                }
                else
                    r = SpatialIndex::Region (reinterpret_cast<const double*>(&i->pos()[0]), reinterpret_cast<const double*>(&i->pos()[0]), this->dimension_);
                this->rtree_->intersectsWithQuery(r, v);
                if (v.intersections()) {
                    i->reload(i->pos() + .02*i->vel(), i->vel());
                }
            }
        });
    }
    centre /= double(this->number_);
    this->setBests_();
    SpatialIndex::Region r;
    if (this->dimension_ == 1){
        double lower[] = {(centre - hyperrect_size)[0], 0};
        double upper[] = {(centre + hyperrect_size)[0], 0};
        r = SpatialIndex::Region(lower, upper, 2);
    }
    else
        r = SpatialIndex::Region(reinterpret_cast<const double*>(&(centre - hyperrect_size)[0]), reinterpret_cast<const double*>(&(centre + hyperrect_size)[0]), this->dimension_);
    IntersectingVisitor v;
    this->rtree_->intersectsWithQuery(r, v);
    if (v.intersections() > this->number_ / 2) {
        this->done_ = true;
    }
}

template <typename T>
T Swarm<T>::bestVal() {
    return this->particles_[0]->getVal();
}

template <typename T>
const std::vector<double>& Swarm<T>::bestX() {
    return this->particles_[0]->pos();
}
