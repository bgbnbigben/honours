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
		mexErrMsgTxt("This function requires 3 arguments\n");
	if (!mxIsNumeric(prhs[0]))
		mexErrMsgTxt("parameter 1 missing!\n");

	// retrieve the heap
	Queue<std::vector<double> >* queue;
	retrieve_queue(prhs[0], queue);

	// delete the queue
	//queue -> ~Queue<std::vector<double> >();
    delete queue;
}
#endif
