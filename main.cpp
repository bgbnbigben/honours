#include <iostream>
#include <iomanip>
//#include <parallel/algorithm>
//#include <parallel/settings.h>
#include <vector>
#include <swarm/swarm.h>
#include <swarm/neldermead.h>
#include <utilities/function.h>
#include <utilities/tsqueue.h>
#include <utilities/vector_ops.h>
#include <utilities/bound.h>
#include <bfgs/bfgs.h>
#include <optionparser.h>
//#include <mpi.h>
constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;

struct Arg : public option::Arg {
    static void printError(const char* msg1, const option::Option& opt, const char* msg2) {
        fprintf(stderr, "%s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
    }

    static option::ArgStatus Unknown(const option::Option& option, bool msg) {
        std::cout << "UNKNOWN CALLED" << std::endl;
        if (msg) printError("Unknown option '", option, "'\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Required(const option::Option& option, bool msg) {
        std::cout << option.arg << " LKSJHDFLKJSHDLFKJ" << std::endl;
        if (option.arg != 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires an argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus NonEmpty(const option::Option& option, bool msg) {
        if (option.arg != 0 && option.arg[0] != 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires a non-empty argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Numeric(const option::Option& option, bool msg) {
        char* endptr = 0;
        if (option.arg != 0 && strtod(option.arg, &endptr)){};
        if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires a numeric argument\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Integer(const option::Option& option, bool msg) {
        char* endptr = 0;
        auto inter = std::bind(strtol, std::placeholders::_1, std::placeholders::_2, 10);
        if (option.arg != 0 && inter(option.arg, &endptr)){};
        if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;

        if (msg) printError("Option '", option, "' requires an integer argument\n");
        return option::ARG_ILLEGAL;
    }
};

enum optionIndex { UNKNOWN, HELP, NUM, DIM, BOUND, NONEMPTY };
const option::Descriptor usage[] = {
    { UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: ./particle [options]\n\n"
                                        "Options: "},
    { HELP, 0, "h", "help", Arg::None,  "  \t-h|--help \tPrint usage and exit."
                                                   },
    { BOUND, 0, "l", "left-window", Arg::Numeric, " \t-l|--left-window "
                             "<double> the left window of the search space."
                             " This will be considered a hard bound." },
    { BOUND, 0, "r", "right-window", Arg::Numeric, " \t-r|--right-window "
                             "<double> the right window of the search space."
                             " This will be considered a hard bound." },
    { NUM, 0, "n", "num", Arg::Integer, " \t-n|--num <int> the number of"
                                        " particles in the swarm." },
    { DIM, 0, "d", "dim", Arg::Integer, " \t-d|--dim <int> the dimension of"
                                        " the problem to be solved." },
    { UNKNOWN, 0, "", "", Arg::None, "\nExamples:\n"
                                             " ./particle -l -10.4 -r 100 "
                                             "-n 10 -d 1024\n"},
    { 0, 0, 0, 0, 0, 0}
};

void invalidInvocation() {
    int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
    option::printUsage(fwrite, stdout, usage, columns);
    exit(0);
}

int main(int argc, char* argv[]) {
    char idstr[32], buff[128];
    char* end = 0;
    int numprocs, rank, i;

    //MPI_Status stat;
    //MPI_Init(&argc, &argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    argv += (argc > 0);
    argc -= (argc > 0);

    option::Stats stats(usage, argc, argv);
    option::Option *buffer = new option::Option[stats.buffer_max];
    option::Option *options = new option::Option[stats.options_max];

    option::Parser parse(usage, argc, argv, options, buffer);
    if (parse.error())
        return 1;

    if (options[HELP] || argc == 0 || options[BOUND].count() != 2 || !options[NUM].count() || !options[DIM].count()) {
        invalidInvocation();
    }

    double left = strtod(options[BOUND].arg, &end);
    double right = strtod(options[BOUND].next()->arg, &end);
    int n = strtol(options[NUM].arg, &end, 10);
    int dim = strtol(options[DIM].arg, &end, 10);;

    if (std::isinf(left)  || std::isnan(left) ||
        std::isinf(right) || std::isnan(right))
        invalidInvocation();

    std::vector<Bound<double>> bounds(n);
    # pragma omp parallel for
        for (int i = 0; i < n; i++) {
            bounds[i].lower = left;
            bounds[i].upper = right;
        }

    if (rank == 0) {
        /* Do master shit. Divide up work, set bounds, etc. */
    } else {
    }
    //__gnu_parallel::_Settings s;
    //s.algorithm_strategy = __gnu_parallel::force_parallel;
    //__gnu_parallel::_Settings::set(s);

    tsqueue<std::vector<double>> q;
    Rosenbrock<double> testFunction;
    BFGS<double, Rosenbrock<double> > bfgs;
        return 0;
    Swarm swarm(&testFunction, n, dim, bounds);

    double bestF = swarm.bestVal();
    std::cout << "Starting at " << swarm.bestX() << std::endl;
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(10);
    while (!swarm.done()) {
        swarm.dance();
        std::cout << "Best in swarm: " << swarm.bestVal() << std::endl;
        std::cout << "Best so far: " << bestF << std::endl;
        if (bestF > swarm.bestVal())
            bestF = swarm.bestVal();
    }
    std::cout << "Swarm took us to " << swarm.bestX() << std::endl;
    q.push(swarm.bestX());

    while (!q.empty()) {
        auto guess = *q.pop();
        std::cout << "Guessing!" << std::endl;
        bfgs.optimize(testFunction, guess, bounds);
        if (bfgs.isSuccess()) {
            std::cout << setiosflags(std::ios::fixed) << std::setprecision(6)
                      << "The optimal value is at " <<  bfgs.getOptValue()
                      << std::endl;
            std::cerr << "f(x) = " << testFunction(bfgs.getOptValue())
                      << std::endl;
        } else {
            std::cout << "************\nBFGS Failed!\n************"
                      << std::endl;
            std::cout << "We started from\n"
                      << swarm.bestX()
                      << std::endl;
        }
    }

    std::cout << "We had " << swarm.numIterations() << " iterations of the swarm"
              << " and " << bfgs.numIterations() << " BFGS iterations"
              << std::endl;

    std::cout << "This is a total of " << testFunction.numCalls()
              << " function calls" << std::endl;
    return 0;
}
