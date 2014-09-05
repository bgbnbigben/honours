#ifndef SWARM_H
#define SWARM_H
#include <vector>
#include <numeric>
#include <utilities/function.h>
#include <utilities/bound.h>
#include <spatialindex/SpatialIndex.h>

class Particle;

class Swarm {
        Function<double>* f_;
        unsigned number_;
        int dimension_;
        SpatialIndex::IStorageManager* memoryStorage_;
        SpatialIndex::StorageManager::IBuffer* memoryBuffer_;
        SpatialIndex::ISpatialIndex* rtree_;
        bool done_;
        int same_;
        unsigned iterations_;
        std::vector<Bound<double>> bounds_;

    public:
        std::vector<Particle*> particles_;
    private:

        void setBests_();

    public:
        Swarm(Function<double>*, int, int, std::vector<Bound<double>>);
        ~Swarm();

        void dance();
        double bestVal();
        const std::vector<double>& bestX();
        bool done() const { return this->done_; }

        int numIterations() const { return this->iterations_; }
        int numParticles() const { return this->number_; }

        static constexpr double c1 = 2.8;
        static constexpr double c2 = 1.3;
        static constexpr double momentum = 0.729844;

};

#endif
