#include "../nomad.3.5.1/src/nomad.hpp"
#include "optionparser.h"
#include "default_eval.h"
#include <cmath>
#include <string>
#include <algorithm>
using namespace std;

const int NUM_VARIABLES = 5;
const int NUM_ITERATIONS = 10000;


struct Arg: public option::Arg
{
  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;

    if (msg) {
        fprintf(stderr, "Option '");
        fwrite(option.name, option.namelen, 1, stderr);
        fprintf(stderr, "' requires a non-empty argument\n");
    }
    return option::ARG_ILLEGAL;
  }
};

enum optionIndex {PROB_TYPE};
const option::Descriptor usage[] = {
    {PROB_TYPE, 0, "p", "problem", Arg::NonEmpty,
        "  -p <name>, \t--problem=<name>. Defaults to Griewangk."},
    {0, 0, 0, 0, 0, 0}
};

class Griewangk : public NOMAD::Evaluator {
    public:
        Griewangk (const NOMAD::Parameters &p) : NOMAD::Evaluator(p) {}

        /* Attempt to optimise
         * f(x) = sum(x^2) / 4000 - prod(cos(x / sqrt(n))) + 1
         */
        bool eval_x (NOMAD::Eval_Point &x, const NOMAD::Double &h_max,
                     bool &count_eval) const {

            const int dimension = this->_p.get_dimension();
            // TODO do this better
            NOMAD::Double prod(1.0);
            for (int i = 0; i < dimension; i++) {
                prod *= cos((x[i] / sqrt(dimension)).value());
            }

            x.set_bb_output(0, x.dot_product(x) / 4000.0 - prod + 1.0);
            count_eval = true;
            return true;
        }
};

int main (int argc, char** argv) {
    option::Stats  stats(usage, argc-1, argv+1);
    option::Option* options = new option::Option[stats.options_max];
    option::Option* buffer  = new option::Option[stats.buffer_max];
    option::Parser parse(usage, argc-1, argv+1, options, buffer);

    NOMAD::Display out(std::cout);
    out.precision(NOMAD::DISPLAY_PRECISION_STD);
    NOMAD::Evaluator *my_eval;

    try {
        NOMAD::begin(argc, argv);
        NOMAD::Parameters p(out);

        p.set_DISPLAY_STATS("bbe ( sol ) obj");
        //p.set_DISPLAY_DEGREE(FULL_DISPLAY);
        p.set_DIMENSION(NUM_VARIABLES);
        p.set_MAX_BB_EVAL(NUM_ITERATIONS);
        // p.set_TMP_DIR("/tmp"); // Defaults to pwd

        vector<NOMAD::bb_output_type> bbot(3);  // definition of
        bbot[0] = NOMAD::OBJ;                   // output types
        bbot[1] = NOMAD::PB;
        bbot[2] = NOMAD::EB;
        p.set_BB_OUTPUT_TYPE(bbot);

        p.set_X0(NOMAD::Point(NUM_VARIABLES, 0.4));  // starting point

        p.set_LOWER_BOUND(NOMAD::Point(NUM_VARIABLES, -6.0)); // all var. >= -6
        NOMAD::Point ub(NUM_VARIABLES);                    // x_4 and x_5 have no bounds
        ub[0] = 5.0;                              // x_1 <= 5
        ub[1] = 6.0;                              // x_2 <= 6
        ub[2] = 7.0;                              // x_3 <= 7
        p.set_UPPER_BOUND(ub);

        p.check();
        
        string data(options[PROB_TYPE].arg?options[PROB_TYPE].arg:"");
        transform(data.begin(), data.end(), data.begin(), ::tolower);
        if (!data.compare("default")) {
            my_eval = new My_Evaluator(p);
        } else {
            my_eval = new Griewangk(p);
        }
        NOMAD::Mads mads(p, my_eval);
        mads.run();
    } catch (exception &e) {
        cerr << "\nNOMAD has been interrupted(" << e.what() << ")\n\n";
    }

    NOMAD::Slave::stop_slaves(out);
    NOMAD::end();

    return EXIT_SUCCESS;
}
