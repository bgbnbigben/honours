/* N.B. This code assumes there are at most as many bounds as particles in the
 * swarm.
 * This assumption is fucked
 */
#include <iostream>
#include <iomanip>
//#include <parallel/algorithm>
//#include <parallel/settings.h>
#include <vector>
#include <limits>
#include <swarm/swarm.h>
#include <swarm/neldermead.h>
#include <utilities/RTreeUtilities.h>
#include <utilities/function.h>
#include <utilities/tsqueue.h>
#include <utilities/vector_ops.h>
#include <utilities/bound.h>
#include <bfgs/bfgs.h>
#include <optionparser.h>
#include <mpi.h>
#include <csignal>
#include <chrono>
#include <string>

#define BOUND_TAG 0x626f756e
#define DIE_TAG 0xd1ed1e
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;
constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;

int current_proc;
int numprocs;

void walltime(int param) {
    std::cerr << "Exceeded walltime with param " << param << std::endl;
}

struct Arg : public option::Arg {
    static void printError(const char* msg1, const option::Option& opt, const char* msg2) {
        fprintf(stderr, "%s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
    }

    static option::ArgStatus Unknown(const option::Option& option, bool msg) {
        if (msg) printError("Unknown option '", option, "'\n");
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus Required(const option::Option& option, bool msg) {
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

    static option::ArgStatus String(const option::Option&, bool) {
        /* By default everything is a string.... */
        return option::ARG_OK;
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

enum optionIndex { UNKNOWN, HELP, NUM, DIM, BOUND, PARTITIONS, NONEMPTY, SIF, NAME };
const option::Descriptor usage[] = {
    { UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: ./particle [options]\n\n"
                                        "Options: "},
    { HELP, 0, "h", "help", Arg::None,  "  \t-h|--help \tPrint usage and exit."
                                                   },
    { NUM, 0, "n", "num", Arg::Integer, " \t-n|--num <int> the number of"
                                        " particles in the swarm." },
    { BOUND, 0, "l", "left-window", Arg::Numeric, " \t-l|--left-window "
                             "<double> the left window of the search space."
                             " This will be considered a hard bound. This"
                             " argument is ignored if using a problem from"
                             " CUTEst." },
    { BOUND, 0, "r", "right-window", Arg::Numeric, " \t-r|--right-window "
                             "<double> the right window of the search space."
                             " This will be considered a hard bound. This"
                             " argument is ignored if using a problem from"
                             " CUTEst." },
    { DIM, 0, "d", "dim", Arg::Integer, " \t-d|--dim <int> the dimension of"
                                        " the problem to be solved. This" 
                                        " argument is ignored if using a"
                                        " problem from CUTEst." },
    { NAME, 0, "a", "name", Arg::String, " \t-a|--name <name> the name of the"
                                         " (C++) function to optimise on."
                                         " Current valid options are "
                                         " 'Paraboloid' or 'Rosenbrock'." },
    { SIF, 0, "s", "sif", Arg::String, " \t-s|--sif <name> the name of the SIF"
                                       " file to test with, from the CUTEst"
                                       " dataset. Your environment must be set"
                                       " up properly for this argument to work."
    },
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

void determinePartitions(int size, int requiredPartitions, long long& actual_, int& n_, int& k_) {
    n_ = std::exp(std::ceil(std::log(requiredPartitions)/size));
    auto num = std::log(requiredPartitions) - size*std::log(n_), denom = std::log(n_ - 1) - std::log(n_);
    k_ = std::floor((double)num/(double)denom);
    actual_ = std::pow(n_ - 1, k_) * std::pow(n_, size - k_);
}

std::vector<std::vector<Bound<double>>> gBounds;

void sendPartitions(std::vector<Bound<double>> bounds, int k, int n, int idx) {
    if (idx == 0) 
        for (unsigned i = 0; i < bounds.size(); i++) {
            bounds[i].variable = i;
            bounds[i].type = Bound<double>::BOTH;
        }
    if (idx == bounds.size()) {
        if (gBounds.size() > 20)
            gBounds.erase(gBounds.begin(), gBounds.begin() + gBounds.size() - 15);
        gBounds.push_back(bounds);
        MPI::Request req = MPI::COMM_WORLD.Isend(&gBounds.back().front(), bounds.size()*sizeof(bounds[0]), MPI::CHAR, (current_proc++ % (numprocs - 1)) + 1, BOUND_TAG);
        return;
    }
    for (int i = 0; i < (n + (idx >= k)); i++) {
        auto newBounds = bounds;
        newBounds[idx].lower = bounds[idx].lower + (i) * (bounds[idx].upper - bounds[idx].lower) / (double)(n + (idx >= k));
        newBounds[idx].upper = bounds[idx].lower + (i + 1)*(bounds[idx].upper - bounds[idx].lower) / (double)(n + (idx >= k));
        sendPartitions(newBounds, k, n, idx+1);
    }
}

int main(int argc, char* argv[]) {
try {
    auto start = std::chrono::high_resolution_clock::now();
    signal(140, walltime);
    MPI::Init(argc, argv);
    char* end = 0;
    int rank;
    double left, right;
    int n, dim;
    long long partitions = 0;
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

    if (options[HELP] || argc == 0 || !options[NUM].count() || ((options[BOUND].count() != 2 || !options[DIM].count()) && options[NAME]) || (!options[SIF] && !options[NAME])) {
        invalidInvocation();
    }

    Function<double>* testFunction;
    BFGS<double> bfgs;

    if (options[SIF]) {
        //testFunction = new CUTEst<double>();
        //dim = (reinterpret_cast<CUTEst<double>*>(testFunction))->getDim();
        //bounds.reserve(dim);
        //bounds.resize(dim);
        //if (rank == 0) {
        //    auto lower = (reinterpret_cast<CUTEst<double>*>(testFunction))->getLowerBounds();
        //    auto upper = (reinterpret_cast<CUTEst<double>*>(testFunction))->getUpperBounds();
        //    for (int i = 0; i < dim; i++) {
        //        std::cout << lower[i] << " <= x_" << i << " <= " << upper[i] << std::endl;
        //        bounds[i].lower = lower[i];
        //        bounds[i].upper = upper[i];
        //        //bounds[i].variable = i;
        //    }
        //}
    } else {
        dim = strtol(options[DIM].arg, &end, 10);
        bounds.reserve(dim);
        bounds.resize(dim);
        if (rank == 0) {
            /* Divide up work, set bounds, etc. */
            left = strtod(options[BOUND].arg, &end);
            right = strtod(options[BOUND].next()->arg, &end);
            if (std::isinf(left)  || std::isnan(left) ||
                std::isinf(right) || std::isnan(right))
                invalidInvocation();
            # pragma omp parallel for
                for (int i = 0; i < dim; i++) {
                    bounds[i].lower = left;
                    bounds[i].upper = right;
                }
        }
        if (std::string(options[NAME].arg) == "Rosenbrock")
            testFunction = new Rosenbrock<double>();
        if (std::string(options[NAME].arg) == "Paraboloid")
            testFunction = new Paraboloid<double>();
    }
    if (!options[PARTITIONS])
        partitions = 2*numprocs;
    else
        partitions = strtol(options[PARTITIONS].arg, &end, 10);

    std::clog << "Expected " << partitions << " partitions to be created" << std::endl;

    n = strtol(options[NUM].arg, &end, 10);
    //n = std::min(n, dim);
    std::clog << "There will be " << n << "particles" << std::endl;

    if (rank == 0) {
        current_proc = 1;
        std::clog << "producing partitions" << std::endl;
        int low, num;
        determinePartitions(bounds.size(), partitions, partitions, low, num);
        std::clog << "We have a total of " << partitions << " partitions and there are " << num << " of size " << low - 1 << std::endl;
        sendPartitions(bounds, num, low - 1, 0);
        std::clog << "Pushed all bounds!" << std::endl;
        assert(bounds.size() == dim);
        gBounds.clear();
        # pragma omp parallel for
            for (int i = 1; i < numprocs; i++) {
                std::clog << "Sending DIE tag to " << i << std::endl;
                MPI::COMM_WORLD.Isend(&bounds.front(), dim*sizeof(bounds[0]), MPI::CHAR, i, DIE_TAG);
            }
        std::clog << "Sent die tags" << std::endl;
        int counter = 0;
        unsigned fCalls = 0, gCalls = 0;
        double globalBest = std::numeric_limits<double>::max();
        std::vector<double> globalSoln(dim);
        std::clog << "Waiting to receive " << partitions << " partitions" << std::endl;
        while (counter < partitions) {
            std::clog << "Waiting for partition number " << counter + 1 << std::endl;
            MPI::Status stat;
            double candidateBest;
            unsigned calls, grads;
            std::vector<double> candidateSoln(dim);
            std::clog << "Prepared to receive a best solution" << std::endl;
            MPI::COMM_WORLD.Recv(&candidateBest, 1, MPI::DOUBLE, MPI::ANY_SOURCE, MPI::ANY_TAG, stat);
            std::clog << "Received the candidate best from processor " << stat.Get_source() << std::endl;
            MPI::COMM_WORLD.Recv(&calls, 1, MPI::UNSIGNED, stat.Get_source(), MPI::ANY_TAG);
            MPI::COMM_WORLD.Recv(&grads, 1, MPI::UNSIGNED, stat.Get_source(), MPI::ANY_TAG);
            MPI::COMM_WORLD.Recv(&candidateSoln.front(), dim, MPI::DOUBLE, stat.Get_source(), MPI::ANY_TAG);
            std::clog << "Received the candidate solution " << candidateBest << " with an associated xval" << std::endl;

            if (candidateBest < globalBest) {
                globalBest = candidateBest;
                globalSoln = candidateSoln;
            }
            fCalls += calls;
            gCalls += grads;
            counter++;
        }

        std::cout << "\n\n\nThe best solution has f value " << globalBest << " and solution\n" << globalSoln << std::endl;
        std::cout << "This is a total of " << fCalls << " function evaluations and " << gCalls << " gradient evaluations (if available)" << std::endl;
    } else {
        MPI::Status status;
        while (true) {
            std::clog << "Processor " << rank << " is blocked" << std::endl;
            MPI::COMM_WORLD.Recv(&bounds.front(), dim*sizeof(Bound<double>), MPI::CHAR, 0, MPI::ANY_TAG, status); 
            std::clog << "Processor " << rank << " has received " << dim*sizeof(Bound<double>) << " bytes worth of bounds" << std::endl;
            if (status.Get_tag() == DIE_TAG) break;

            tsqueue<std::vector<double>> q;
            Swarm swarm(testFunction, n, dim, bounds);

            double bestF = swarm.bestVal();
            std::clog.unsetf ( std::ios::floatfield );
            std::clog.precision(10);
            std::clog << "The bounds are: ";
            for (unsigned i = 0; i < dim; i++)
                std::clog << bounds[i].lower << " <= x_" << i << " <= " << bounds[i].upper << std::endl;
            while (!swarm.done()) {
                swarm.dance();
                if (bestF > swarm.bestVal())
                    bestF = swarm.bestVal();
            }
            std::clog << "On processor " << rank << " the swarm took us to " << swarm.bestX() << std::endl;

            q.push(swarm.bestX());

            while (!q.empty()) {
                auto guess = *q.pop();
                // TODO mesh search up in this bitch
                bfgs.optimize(*testFunction, guess, bounds);
                if (bfgs.isSuccess()) {
                    std::clog << setiosflags(std::ios::fixed) << std::setprecision(6)
                              << "The optimal value is at " <<  bfgs.getOptValue()
                              << std::endl;
                    std::clog << "f(x) = " << testFunction->operator()(bfgs.getOptValue())
                              << std::endl;
                } else {
                    std::clog << "************\nBFGS Failed!\n************"
                              << std::endl;
                    std::clog << "We started from\n"
                              << swarm.bestX()
                              << std::endl;
                }
            }

            std::clog << "We had " << swarm.numIterations() << " iterations of the swarm"
                      << " and " << bfgs.numIterations() << " BFGS iterations"
                      << std::endl;

            std::clog << "This is a total of " << testFunction->numCalls()
                      << " function calls" << std::endl;

            auto val = testFunction->operator()(bfgs.getOptValue());
            auto calls = testFunction->numCalls(), grads = testFunction->numGrads();
            MPI::COMM_WORLD.Isend(&val, 1, MPI::DOUBLE, 0, 0);
            MPI::COMM_WORLD.Isend(&calls, 1, MPI::UNSIGNED, 0, 0);
            MPI::COMM_WORLD.Isend(&grads, 1, MPI::UNSIGNED, 0, 0);
            MPI::COMM_WORLD.Isend(&bfgs.getOptValue().front(), dim, MPI::DOUBLE, 0, DIE_TAG);
        }
    }
    
    //__gnu_parallel::_Settings s;
    //s.algorithm_strategy = __gnu_parallel::force_parallel;
    //__gnu_parallel::_Settings::set(s);

    MPI::Finalize();
    std::cout << std::flush;
    if (rank == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Execution took " << std::chrono::duration_cast<std::chrono::hours>(end - start).count() << "h"
                                       << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << "m"
                                       << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s"
                                       << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
                                       << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us"
                                       << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "ns";
    }
    return 0;
} catch (Tools::IllegalArgumentException& e) {
    std::cerr << e.what() << std::endl;
    return -1;
}
}
