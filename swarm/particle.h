#ifndef PARTICLE_H
#define PARTICLE_H
#include <vector>

class Particle {
        std::vector<double> position_;
        std::vector<double> velocity_;
        std::vector<double> bestPosition_;
        double leftWindow_;
        double rightWindow_;

        double cost_;

    public:
        Particle(double, double, int);
        void step(const std::vector<double>&, double, double, double);
        double& operator[](int i);
        inline void setVal(double cost) { 
            if (cost < this->cost_) 
                this->bestPosition_ = this->position_;
            this->cost_ = cost;
        }
        inline double getVal() const { return this->cost_; }
        inline const std::vector<double>& pos() { return this->position_; }
        inline const std::vector<double>& vel() { return this->velocity_; }
        static constexpr double mutation_prob = 0.05;
        void reload(const std::vector<double>&, const std::vector<double>&);
};

#endif
