#include <iostream>
#include <iomanip>
//#include <parallel/algorithm>
//#include <parallel/settings.h>
#include <vector>
#include <limits>
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
#define DIE_TAG 0xd1ed1e
constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;

MPI::Datatype BoundType;
std::vector<std::vector<Bound<double>>> globalPartitions;
int current_proc;
int numprocs;


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

void producePartitions(std::vector<Bound<double>> bounds, unsigned depth, unsigned requiredPartitions) {
    if (requiredPartitions == 1) {
        /* TODO wait for this to buffer */
        globalPartitions.push_back(bounds);
        std::cout << "Pushing vector to " << ((current_proc + 1) % (numprocs - 1)) + 1 << std::endl;
        std::cout << "Sending " << sizeof(bounds[0])*bounds.size() << " of size " << sizeof(MPI::CHAR) << std::endl;
        MPI::Request req = MPI::COMM_WORLD.Isend(&bounds.front(), bounds.size()*sizeof(bounds[0]), MPI::CHAR, (current_proc++ % (numprocs - 1)) + 1, BOUND_TAG);
        return;
    }
    assert(depth <= bounds.size());
    int n; for (n = 2; n < 10000; n++) if (requiredPartitions - pow(std::pow(n, depth), (bounds.size() - (depth - 1))) <= 0) break;
    int extraPartitions = requiredPartitions - n * (requiredPartitions / n);
    requiredPartitions /= (n+extraPartitions);
    for (int i = 0; i < n + extraPartitions; i++) {
        auto newBounds = bounds;
        newBounds[depth-1].lower = bounds[depth-1].lower + (i) * (bounds[depth-1].upper - bounds[depth-1].lower) / (double)(n + extraPartitions);
        newBounds[depth-1].upper = bounds[depth-1].lower + (i + 1)*(bounds[depth-1].upper - bounds[depth-1].lower) / (double)(n + extraPartitions);
        producePartitions(newBounds, depth+1, requiredPartitions);
    }
}

int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    char* end = 0;
    int rank;
    double left, right;
    int n, dim;
    int partitions = 0;
    std::vector<Bound<double>> bounds;

    MPI::Status stat;
    numprocs = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_ARE_FATAL);


    argv += (argc > 0);
    argc -= (argc > 0);

    option::Stats stats(usage, argc, argv);
    option::Option *buffer = new option::Option[stats.buffer_max];
    option::Option *options = new option::Option[stats.options_max];

    option::Parser parse(usage, argc, argv, options, buffer);
    if (parse.error())
        return 1;

    n = strtol(options[NUM].arg, &end, 10);
    dim = strtol(options[DIM].arg, &end, 10);

    if (options[HELP] || argc == 0 || options[BOUND].count() != 2 || !options[NUM].count() || !options[DIM].count()) {
        invalidInvocation();
    }

    if (rank == 0) {
        /* Divide up work, set bounds, etc. */
        left = strtod(options[BOUND].arg, &end);
        right = strtod(options[BOUND].next()->arg, &end);
        if (!options[PARTITIONS])
            partitions = 2*numprocs;
        else
            partitions = strtol(options[PARTITIONS].arg, &end, 10);

        if (std::isinf(left)  || std::isnan(left) ||
            std::isinf(right) || std::isnan(right))
            invalidInvocation();

    }
    bounds.reserve(n);
    bounds.resize(n);

    // http://stackoverflow.com/a/16747896/699674
    const int    nItems=4;
    int          blocklengths[nItems] = {1, 1, 1, 1};
    MPI::Datatype types[nItems] = {MPI::DOUBLE, MPI::DOUBLE, MPI::LONG, MPI::INT};
    MPI::Aint     offsets[nItems];

    offsets[0] = offsetof(Bound<double>, lower);
    offsets[1] = offsetof(Bound<double>, upper);
    offsets[2] = offsetof(Bound<double>, type);
    offsets[3] = offsetof(Bound<double>, variable);

    MPI::Datatype BoundTypeProto = MPI::Datatype::Create_struct(nItems, blocklengths, offsets, types);

    // Get the constructed type lower bound and extent
    MPI::Aint lb, extent;
    BoundTypeProto.Get_extent(lb, extent);

    // Get the actual distance between to vector elements
    // (this might not be the best way to do it - if so, substitute a better one)
    extent = (char*)&bounds[1] - (char*)&bounds[0];

    // Create a resized type whose extent matches the actual distance
    BoundType = BoundTypeProto.Create_resized(lb, extent);
    BoundType.Commit();

    if (rank == 0) {
        current_proc = 1;
        # pragma omp parallel for
            for (int i = 0; i < n; i++) {
                bounds[i].lower = left;
                bounds[i].upper = right;
            }
        std::cout << "Attempting to send partitions " << std::endl;
        producePartitions(bounds, 1, partitions);
        std::cout << "Partitions sent!" << std::endl;
        # pragma omp parallel for
            for (int i = 1; i < numprocs; i++) {
                MPI::COMM_WORLD.Send(&bounds.front(), n*sizeof(bounds[0]), MPI::CHAR, i, DIE_TAG);
            }
        std::cout << "Told all the processors to die" << std::endl;
        int counter = 0;
        double globalBest = std::numeric_limits<double>::max();
        std::vector<double> globalSoln(dim);
        while (counter < numprocs - 1) {
            MPI::Status stat;
            double candidateBest;
            std::vector<double> candidateSoln(dim);
            std::cout << "Prepared to receive a best solution" << std::endl;
            MPI::COMM_WORLD.Recv(&candidateBest, 1, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
            std::cout << "Recieved the candidate best" << std::endl;
            MPI::COMM_WORLD.Recv(&candidateSoln.front(), dim, MPI::DOUBLE, stat.Get_source(), MPI::ANY_TAG);
            std::cout << "Received the candidate solution " << candidateBest << " with an associated xval" << std::endl;

            if (candidateBest < globalBest) {
                globalBest = candidateBest;
                globalSoln = candidateSoln;
            }
            counter++;
        }

        std::cout << "\n\n\nThe best solution has f value " << globalBest << " and solution\n" << globalSoln << std::endl;
    } else {
        MPI::Status status;
        while (true) {
            MPI::COMM_WORLD.Recv(&bounds.front(), n*sizeof(Bound<double>), MPI::CHAR, 0, MPI::ANY_TAG, status); 
            std::cout << "Processor " << rank << " has received " << n*sizeof(Bound<double>) << " bytes worth of bounds" << std::endl;
            if (status.Get_tag() == DIE_TAG) break;

            tsqueue<std::vector<double>> q;
            Rosenbrock<double> testFunction;
            BFGS<double, Rosenbrock<double> > bfgs;
            Swarm swarm(&testFunction, n, dim, bounds);

            double bestF = swarm.bestVal();
            std::cout << "Starting at " << swarm.bestX() << std::endl;
            std::cout.unsetf ( std::ios::floatfield );
            std::cout.precision(10);
            while (!swarm.done()) {
                swarm.dance();
                if (bestF > swarm.bestVal())
                    bestF = swarm.bestVal();
                std::cout << bestF << " ";
            }
            std::cout << "On processor " << rank << " the swarm took us to " << swarm.bestX() << std::endl;
            q.push(swarm.bestX());

            while (!q.empty()) {
                auto guess = *q.pop();
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

            auto val = testFunction(bfgs.getOptValue());
            MPI::COMM_WORLD.Isend(&val, 1, MPI::DOUBLE, 0, 0);
            MPI::COMM_WORLD.Isend(&bfgs.getOptValue().front(), dim, MPI::DOUBLE, 0, DIE_TAG);
        }
    }
    
    //__gnu_parallel::_Settings s;
    //s.algorithm_strategy = __gnu_parallel::force_parallel;
    //__gnu_parallel::_Settings::set(s);

    //Free up the type
    BoundType.Free();

    MPI::Finalize();
    return 0;
}
