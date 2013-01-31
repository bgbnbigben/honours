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
	if (nrhs != 1)
		mexErrMsgTxt("This function requires 3 arguments\n");
	if (!mxIsNumeric(prhs[0]))
		mexErrMsgTxt("parameter 1 missing!\n");

	// retrieve the heap
	Queue<std::vector<double> >* queue;
	retrieve_queue(prhs[0], queue);

	// query top element in the PQ
    std::vector<double> curr = queue->peek();

	// return its value in the output
	plhs[0] = mxCreateDoubleMatrix(1, curr.size(), mxREAL);
    for (int i = 0; i < curr.size(); i++) {
        mxGetPr(plhs[0])[i] = curr[i];
    }
}
#endif
