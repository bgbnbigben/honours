#include "Queue.h"
#include <iostream>

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


void retrieve_val(const mxArray* matptr, std::vector<double>& val){
	// check that I actually received something
	if (matptr == NULL)
		mexErrMsgTxt("missing third parameter (element index)\n");

	// retrieve index
    int num_elements = mxGetN(matptr) * mxGetM(matptr);
    val.resize(num_elements);
    for (int i = 0; i < num_elements; i++) {
        val[i] = mxGetPr(matptr)[i];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if (nrhs != 2)
		mexErrMsgTxt("This function requires 2 arguments\n");
	if (!mxIsNumeric(prhs[0]))
		mexErrMsgTxt("parameter 1 missing!\n");
	if (!mxIsNumeric(prhs[1]))
		mexErrMsgTxt("parameter 2 missing!\n");
    
	// retrieve the heap
	Queue<std::vector<double> >*  queue;
	retrieve_queue(prhs[0], queue);

    std::vector<double> val;
	retrieve_val(prhs[1], val);

    queue->push(val);
}
#endif
