#include <iostream>
#include <iomanip>
//#include <parallel/algorithm>
//#include <parallel/settings.h>
#include <vector>
#include <swarm/swarm.h>
#include <swarm/neldermead.h>
#include <utilities/function.h>
#include <utilities/tsqueue.h>
#include <bfgs/bfgs.h>
constexpr double NelderMead::alpha;
constexpr double NelderMead::gamma;
constexpr double NelderMead::rho;
constexpr double NelderMead::sigma;

int main() {
    //NelderMead n(8, 0, 320);
    //n.addSimplexPoint({10, 10, 10, 10, 10, 10, 10, 10});
    //n.addSimplexPoint({-10, -10, -10, -10, -10, -10, -10, -10});
    //n.addSimplexPoint({20, 20, 20, 20, 20, 20, 20, 20});
    //n.addSimplexPoint({-20, -20, -20, -20, -20, -20, -20, -20});
    //n.addSimplexPoint({30, 30, 30, 30, 30, 30, 30, 30});
    //n.addSimplexPoint({-30, -30, -30, -30, -30, -30, -30, -30});
    //n.addSimplexPoint({40, 40, 40, 40, 40, 40, 40, 40});
    //n.addSimplexPoint({-40, -40, -40, -40, -40, -40, -40, -40});
    //n.addSimplexPoint({50, 50, 50, 50, 50, 50, 50, 50});
    //std::vector<double> soln = n.drive(new Function());
    //std::cout << (*(new Function()))(soln) << ":\t";
    //std::for_each(soln.begin(), soln.end(), [](double i) {
    //    std::cout << i << "\t";
    //});
    //std::cout << std::endl;


    //__gnu_parallel::_Settings s;
    //s.algorithm_strategy = __gnu_parallel::force_parallel;
    //__gnu_parallel::_Settings::set(s);


    tsqueue<std::vector<double>> q;
    Rosenbrock<double> testFunction;

    BFGS<double, Rosenbrock<double> > bfgs;
    Swarm swarm(&testFunction, 300, 10);
    double bestF = swarm.bestVal();
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(10);
    while (!swarm.done()) {
        swarm.dance();
        std::cout << "Best in swarm: " << swarm.bestVal() << std::endl;
        std::cout << "Best so far: " << bestF << std::endl;
        if (bestF > swarm.bestVal())
            bestF = swarm.bestVal();
    }
    q.push(swarm.bestX());

    while (!q.empty()) {
        auto guess = *q.pop();
        std::cout << "Guessing!" << std::endl;
        bfgs.optimize(testFunction, guess, 1.0e-10, 1000);
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

    //std::cout << "There are " << swarm.numBlockedOff() << " items in the tabu list" << std::endl;
    return 0;
}
