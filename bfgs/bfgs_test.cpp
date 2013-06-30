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
    std::vector<double> x0(10);
    for (int i = 0; i < 10; i++) x0[i] = 0.0;

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
