#include "Queue.h"
#include <vector>

//------------------------------- MATLAB -------------------------------------//
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define toSysout(...) printf(__VA_ARGS__)
 #define exit_with_error(...)           \
 do {                                   \
		 fprintf(stdout, "Error: ");    \
		 fprintf(stdout, __VA_ARGS__ ); \
		 fprintf(stdout, "\n" );        \
		 exit(1);                       \
 } while(0)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if (nrhs != 1)
		mexErrMsgTxt("This function requires an argument\n");
	if (!mxIsNumeric(prhs[0]))
		mexErrMsgTxt("parameter 1 missing!\n");

	// retrieve the queue
	Queue<std::vector<double> >*  queue;
	retrieve_queue(prhs[0], queue);

	std::vector<double> curr = queue->pop();
	plhs[0] = mxCreateDoubleMatrix(1, curr.size(), mxREAL);
    for (int i = 0; i < curr.size(); i++) {
        mxGetPr(plhs[0])[i] = curr[i];
    }
	//*mxGetPr(plhs[0]) = curr;

}
#endif
