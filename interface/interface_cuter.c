/* ====================================================
 * CUTEst interface simulating a black box for NOMAD.
 * April 25, 2013
 *
 * D. Orban from an earlier version by S. Le Digabel.
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#ifdef FLT_MAX
#define INFINITY FLT_MAX
#else
#define INFINITY 1.0e+20
#endif

#include "cutest.h"

int CUTEst_nvar;        /* number of variables */
int CUTEst_ncon;        /* number of constraints */

int main(int argc , char** argv) {
    char* fname = "OUTSDIF.d"; /* CUTEst data file */
    int funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    int iout  = 6;         /* FORTRAN unit number for error output */
    int io_buffer = 11;    /* Internal input/output buffer */
    int ierr;              /* Exit flag for various calls */

    double *x, *x_l, *x_u; /* position, lower bound, upper bound */
    double *grad;
    double *y = NULL, *y_l = NULL, *y_u = NULL;
    bool *equatn = NULL, *linear = NULL;
    int eqn_order = 0, l_order = 0, v_order = 0;
    bool constrained = false;

    double obj;
    int i;

    /* Open problem description file OUTSDIF.d */
    ierr = 0;
    FORTRAN_open(&funit, fname, &ierr);
    if (ierr) {
        return -1;
    }

    /* Determine problem size */
    CUTEST_cdimen(&ierr, &funit, &CUTEst_nvar, &CUTEst_ncon);
    if (ierr) {
        return -2;
    }

    /* Determine whether to call constrained or unconstrained tools */
    if (CUTEst_ncon) constrained = true;

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEst */
    x = (double*)malloc((CUTEst_nvar)*sizeof(double));
    x_l = (double*)malloc((CUTEst_nvar)*sizeof(double));
    x_u = (double*)malloc((CUTEst_nvar)*sizeof(double));
    grad = (double*)malloc((CUTEst_nvar)*sizeof(double));

    if (constrained) {
        equatn = (bool*)malloc((CUTEst_ncon+1)*sizeof(bool));
        linear = (bool*)malloc((CUTEst_ncon+1)*sizeof(bool));
        y = (double*)malloc((CUTEst_ncon+1)*sizeof(double));
        y_l = (double*)malloc((CUTEst_ncon+1)*sizeof(double));
        y_u = (double*)malloc((CUTEst_ncon+1)*sizeof(double));
        CUTEST_csetup(&ierr, &funit, &iout, &io_buffer,
                &CUTEst_nvar, &CUTEst_ncon,
                x, x_l, x_u, y, y_l, y_u,
                equatn, linear, &eqn_order, &l_order, &v_order);
        if (ierr) {
            return -3;
        }
    } else {
        equatn = (bool*)malloc((1)*sizeof(bool));
        linear = (bool*)malloc((1)*sizeof(bool));
        y_l = (double*)malloc((1)*sizeof(double));
        y_u = (double*)malloc((1)*sizeof(double));
        CUTEST_usetup(&ierr, &funit, &iout, &io_buffer,
                &CUTEst_nvar, x, x_l, x_u);
        if (ierr) {
            return -3;
        }
    }

    FORTRAN_close(&funit, &ierr);

    /* If an input vector is supplied, use it.
     * Otherwise, use the problem's initial guess. */
    if (argc >= 2) {
        /* See if problem dimension is requested */
        if (strcmp(argv[1], "--nvar") == 0) {
            printf("%d\n", CUTEst_nvar);
            return 0;
        }

        for (i = 0 ; i < CUTEst_nvar; i++)
            scanf("%lf" , &x[i]);
    }

    printf("Starting from the point:\n");
    for (i = 0; i < CUTEst_nvar; i++)
        printf("%lf\n", x[i]);
                
    printf("With bounds:\n");
    for (i = 0; i < CUTEst_nvar; i++)
        printf("%lf <= x_%d <= %lf\n", x_l[i], i, x_u[i]);

    if (constrained) {
        /* Recycle the array y to store constraint values */
        //CUTEST_cfn(&ierr, &CUTEst_nvar, &CUTEst_ncon, x, &obj, y);
        logical a = true;
        CUTEST_cofg(&ierr, &CUTEst_nvar, x, &obj, grad, &a);
        if (ierr) {
            return -4;
        }
        printf("%21.15e ", obj);
        for (i = 0 ; i < CUTEst_ncon ; i++)
            printf("%21.15e ", y[i]);
        printf("\n");
        /* n, m, X, Y, gradient_of_lagrangian, grad, j_trans, l_1, l_2, J */
        /* J = double[l_1][l_2], if j_trans is true this is transposed, and J
         * is the Jacobian matrix of constraint functions
         */
        //CUTEST_cgr(&ierr, &CUTEst_nvar, &CUTEst_ncon, x, y, 
        //           &false, grad, &false, n, m, jac);
    } else {
        logical a = true;
        CUTEST_uofg(&ierr, &CUTEst_nvar, x, &obj, grad, &a);
        if (ierr) {
            return -4;
        }
        printf("%21.15e\n", obj);
        printf("Gradient = \n");
        for (i = 0; i < CUTEst_nvar; i++)
            printf("%lf\n", grad[i]);
    }

    /* Free workspace */
    free(x);
    free(x_l);
    free(x_u);
    free(grad);
    free(y);
    free(y_l);
    free(y_u);
    free(equatn);
    free(linear);

    return 0;
}

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
