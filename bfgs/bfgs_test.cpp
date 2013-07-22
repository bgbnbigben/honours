/*****************************************************************************
 *                               bfgs_test.cpp
 *
 * BFGS method testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <vector>
#include "bfgs.h"
#include <utilities/function.h>

int main() {
    Rosenbrock<double> f;
    //std::vector<double> x0 = {-0.9734268716, 0.9629549863, 0.944658782,   0.8867037531,   0.7815048012,
    //                           0.600754717,  0.3573148646, 0.1168369113, -0.001853005818, 0.05175628025};
    std::vector<double> x0(10);

    BFGS< double, Rosenbrock<double> > bfgs;
    bfgs.optimize( f, x0, double(1.0e-10), 1000);
    std::vector<double> xmin = bfgs.getOptValue();
    if (bfgs.isSuccess()) {
        int N = bfgs.numIterations();
        std::cout << "The iterative number is:  " << N << "\n";
        std::cout << "\nThe number of function calculation is:   "
                  << f.numCalls() << "\n";
        std::cout << "\nThe gradient's norm at x is:   "
                  << bfgs.getGradNorm()[N] << "\n";
    } else {
        std::cout << "The optimal solution can't be found!";
    }
    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    std::cout << "\nThe value of x is:   " << xmin;
    std::cout << "\nThe minimum value of f(x) is:   " << f(xmin) << "\n" << std::endl;


    return 0;
}
