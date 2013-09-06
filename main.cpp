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
constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;

int main() {
    //__gnu_parallel::_Settings s;
    //s.algorithm_strategy = __gnu_parallel::force_parallel;
    //__gnu_parallel::_Settings::set(s);

    tsqueue<std::vector<double>> q;
    Rosenbrock<double> testFunction;
    BFGS<double, Rosenbrock<double> > bfgs;
    Swarm swarm(&testFunction, 10, 10);

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
        bfgs.optimize(testFunction, guess, 1.0e-10, 5000);
        if (bfgs.isSuccess()) {
            std::cout << setiosflags(std::ios::fixed) << std::setprecision(6)
                      << "The optimal value is at " <<  bfgs.getOptValue()
                      << std::endl;
            std::cout << "f(x) = " << testFunction(bfgs.getOptValue())
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
