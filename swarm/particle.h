#ifndef PARTICLE_H
#define PARTICLE_H
#include <vector>
//#include <utilities/bound.h>
#include <utilities/Variable.h>

template <typename T>
class Particle {
        std::vector<Variable> velocity_;
        std::vector<Variable> bestPosition_;
        T cost_;

        VariableType getRandomInWindow_(int i);

    public:
        // TODO wrap these into the one.
        std::vector<VariableContainer> position_;
        Particle(const std::vector<VariableContainer>&, int);
        void step(const std::vector<Variable>&, double, double, double);
        void clamp();
        Variable operator[](int i);
        inline void setVal(T cost) { 
            if (cost < this->cost_) 
                this->bestPosition_ = this->position_;
            this->cost_ = cost;
        }
        inline T getVal() const { return this->cost_; }
        inline const std::vector<Variable>& pos() { return this->position_; }
        inline const std::vector<Variable>& vel() { return this->velocity_; }
        static constexpr double mutation_prob = 0.05;
        void reload(const std::vector<Variable>&, const std::vector<Variable>&);
};

#endif
