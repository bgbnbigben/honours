#include "../nomad.3.5.1/src/nomad.hpp"
class My_Evaluator : public NOMAD::Evaluator {
    public:
        My_Evaluator (const NOMAD::Parameters &p) : NOMAD::Evaluator(p) {}

        /* Claims to attempt to optimise
         * f(x) = x_5, such that 
         * -25 + sum_{i=1}^{5} (x_i - 1)^2 \leq 0
         *  25 - sum_{i=1}^{5} (x_i + 1)^2 \leq 0
         * x_i \geq -6, x_1 \leq 5, x_2 \leq 6, x_3 \leq 7
         */
        bool eval_x (NOMAD::Eval_Point&, const NOMAD::Double&, bool&) const;
};

