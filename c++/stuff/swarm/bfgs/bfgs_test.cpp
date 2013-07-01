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
using namespace std;


int main() {
    Rosenbrock<double> f;
    std::vector<double> x0 = {-0.9734268716,
        0.9629549863,
        0.944658782,
        0.8867037531,
        0.7815048012,
        0.600754717,
        0.3573148646,
        0.1168369113,
        -0.001853005818,
        0.05175628025};

    BFGS< double, Rosenbrock<double> > bfgs;
    bfgs.optimize( f, x0, double(1.0e-10), 1000);
    if( bfgs.isSuccess() )
    {
        std::vector<double> xmin = bfgs.getOptValue();
        int N = bfgs.getItrNum();
        cout << "The iterative number is:   " << N << endl << endl;
        cout << "The number of function calculation is:   "
             << bfgs.getFuncNum() << endl << endl;
        cout << setiosflags(ios::fixed) << setprecision(4);
        cout << "The optimal value of x is:   " << xmin << endl;
        cout << "The minimum value of f(x) is:   " << f(xmin) << endl << endl;
        cout << "The gradient's norm at x is:   "
             << bfgs.getGradNorm()[N] << endl << endl;
    }
    else
        cout << "The optimal solution  cann't be found!" << endl;

    return 0;
}
