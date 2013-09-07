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
#include <mpi.h>

#define BOUND_TAG 0x626f756e
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

    static option::ArgStatus OptionalInteger(const option::Option& option, bool msg) {
        return Arg::Integer(option, msg);
    }
};

enum optionIndex { UNKNOWN, HELP, NUM, DIM, BOUND, PARTITIONS, NONEMPTY };
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
    { PARTITIONS, 0, "p", "partitions", Arg::OptionalInteger,
        " \t-p|--partitions <int> The number of partitions to split the search"
        " space into. The more partitions, the higher the probability of "
        " finding the global optimum, however this will slow the search down."
        " The default value is twice the available number of MPI nodes." },
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
    char* end = 0;
    int numprocs, rank;
    double left, right;
    int n, dim;
    int partitions = 0;
    std::vector<Bound<double>> bounds;

    MPI::Status stat;
    MPI::Init(argc, argv);
    numprocs = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();

    if (rank == 0) {
        /* Divide up work, set bounds, etc. */
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

        left = strtod(options[BOUND].arg, &end);
        right = strtod(options[BOUND].next()->arg, &end);
        n = strtol(options[NUM].arg, &end, 10);
        dim = strtol(options[DIM].arg, &end, 10);

        if (!options[PARTITIONS])
            partitions = 2*numprocs;
        else
            partitions = strtol(options[PARTITIONS].arg, &end, 10);

        if (std::isinf(left)  || std::isnan(left) ||
            std::isinf(right) || std::isnan(right))
            invalidInvocation();

    }
    MPI::COMM_WORLD.Bcast(&n, 1, MPI_INT, 0);
    MPI::COMM_WORLD.Bcast(&dim, 1, MPI_INT, 0);
    bounds.reserve(n);
    bounds.resize(n);

    // http://stackoverflow.com/a/16747896/699674
    const int    nItems=2;
    int          blocklengths[nItems] = {1, 1};
    MPI::Datatype types[nItems] = {MPI::DOUBLE, MPI::DOUBLE};
    MPI::Datatype MPI_BoundType_proto, MPI_BoundType;
    MPI::Aint     offsets[nItems];

    offsets[0] = offsetof(Bound<double>, lower);
    offsets[1] = offsetof(Bound<double>, upper);

    MPI::Datatype BoundTypeProto = MPI::Datatype::Create_struct(nItems, blocklengths, offsets, types);

    // Get the constructed type lower bound and extent
    MPI::Aint lb, extent;
    BoundTypeProto.Get_extent(lb, extent);

    // Get the actual distance between to vector elements
    // (this might not be the best way to do it - if so, substitute a better one)
    extent = (char*)&bounds[1] - (char*)&bounds[0];

    // Create a resized type whose extent matches the actual distance
    MPI::Datatype BoundType = BoundTypeProto.Create_resized(lb, extent);
    BoundType.Commit();

    if (rank == 0) {
        # pragma omp parallel for
            for (int i = 0; i < n; i++) {
                bounds[i].lower = left;
                bounds[i].upper = right;
            }
        int splits = partitions/n;
        int remainder = partitions - splits*n;

        int current_proc = 1;
        for (int i = 0; i < n; i++, remainder--) {
            // need to divide bound[i] up into splits+remainder partitions.
            double lower = bounds[i].lower, upper = bounds[i].upper;
            for (int j = 0; j < splits + (partitions > 0); j++) {
                if (j == 0) {
                    bounds[i].lower = lower;
                } else {
                    bounds[i].lower = lower + (j)*(upper-lower) / (double)(splits+(partitions>0));
                }
                if (j == splits - 1) {
                    bounds[i].upper = upper;
                } else {
                    bounds[i].upper = lower + (j+1)*(upper-lower) / (double)(splits+(partitions>0));
                }

                MPI::COMM_WORLD.Send(&bounds.front(), n, BoundType, current_proc++, BOUND_TAG);
            }
        }
    } else {
        MPI::COMM_WORLD.Recv(&bounds.front(), n, BoundType, 0, BOUND_TAG);
    }
    
    //Free up the type
    BoundType.Free();

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
    MPI::Finalize();
    return 0;
}
