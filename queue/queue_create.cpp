#include "Queue.h"
#include <vector>

//------------------------------- MATLAB -------------------------------------//
 #define toSysout(...) printf(__VA_ARGS__)
 #define exit_with_error(...)           \
 do {                                   \
		 fprintf(stdout, "Error: ");    \
		 fprintf(stdout, __VA_ARGS__ ); \
		 fprintf(stdout, "\n" );        \
		 exit(1);                       \
 } while(0)
#ifdef MATLAB_MEX_FILE
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// instantiate the queue
	Queue<std::vector<double> >* q = new Queue<std::vector<double> >();

	// convert the points to double
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* pointer_to_queue = mxGetPr(plhs[0]);
	pointer_to_queue[0] = (long) q;
}
#endif
