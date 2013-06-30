#include <iostream>
#include <iomanip>
//#include <parallel/algorithm>
//#include <parallel/settings.h>
#include <vector>
#include <swarm/swarm.h>
#include <swarm/neldermead.h>
#include <utilities/function.h>
#include <utilities/tsqueue.h>
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

    Swarm swarm(&testFunction, 300, 10);
    double bestF = swarm.bestVal();
    int same = 0;
    std::cout.unsetf ( std::ios::floatfield );
    std::cout.precision(10);
    while (same < 100) {
        swarm.dance();
        std::cout << "Best in swarm: " << swarm.bestVal() << std::endl;
        std::cout << "Best so far: " << bestF << std::endl;
        if (bestF > swarm.bestVal()) {
            bestF = swarm.bestVal();
            same = 0;
        } else {
            same++;
        }
    }
    q.push(swarm.bestX());

    while (!q.empty()) {
        auto guess = q.pop();
    }

    std::cout << bestF << ":\t( ";
    std::for_each(swarm.bestX().begin(), swarm.bestX().end(), [](double i) {
        std::cout << i << ",\t";
    });
    std::cout << ")\n" << std::endl;

    //std::cout << "There are " << swarm.numBlockedOff() << " items in the tabu list" << std::endl;
    return 0;
}
