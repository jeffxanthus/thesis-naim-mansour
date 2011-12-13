/*
 * TOMLAB SA SHARED LIBRARY EXAMPLE
 * application.c
 *
 * If the preprocessor flag MCC4 is defined, then compile code
 * supported by MCC4, else compile code supported by MCC3.
 *
 *
 * Fredrik Hellman, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
 * Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.5.0$
 * Written Nov 11, 2004.   Last modified Nov 11, 2004.
 */

#include "libcplexqp.h"
#include <stdio.h>

int main()
{
  mxArray *F, *c, *A, *b_L, *b_U, *x_L, *x_U;
  mxArray *x, *f;

  double Fdata[2*2] = {2, 2, 2, 2};
  double cdata[2] = {2, 6};
  double Adata[3*2] = {1, 1, 0, 1, 0, 1};
  double b_Ldata[3] = {3, 0, 0};
  double b_Udata[3] = {1e20, 1e20, 1e20};
  double x_Ldata[2] = {-1, -1};
  double x_Udata[2] = {100, 100};

#ifdef MCC4
  printf("Initializing application...\n");

  /* This initialization must be made at the beginning in order
   * to be able to use the mx/mex routines */

  mclInitializeApplication(NULL,0);
#endif

  printf("Initializing stand alone library...\n");

#ifdef MCC4
  if (!libcplexqpInitialize())
  {
    fprintf(stderr,"Could not initialize the cplexqp library.");
    return -1;
  } 
#else
  libcplexqpInitialize();
#endif

  printf("Creating user data matrices...\n");

  /* Create MATLAB input matrices */
  F   = mxCreateDoubleMatrix(2, 2, mxREAL);
  c   = mxCreateDoubleMatrix(2, 1, mxREAL);
  A   = mxCreateDoubleMatrix(3, 2, mxREAL);
  b_L = mxCreateDoubleMatrix(3, 1, mxREAL);
  b_U = mxCreateDoubleMatrix(3, 1, mxREAL);
  x_L = mxCreateDoubleMatrix(2, 1, mxREAL);
  x_U = mxCreateDoubleMatrix(2, 1, mxREAL);

  /* Copy the user data into the MATLAB matrices */
  memcpy(mxGetPr(F)  ,   Fdata, 2*2*sizeof(double));
  memcpy(mxGetPr(c)  ,   cdata, 2*1*sizeof(double));
  memcpy(mxGetPr(A)  ,   Adata, 3*2*sizeof(double));
  memcpy(mxGetPr(b_L), b_Ldata, 3*1*sizeof(double));
  memcpy(mxGetPr(b_U), b_Udata, 3*1*sizeof(double));
  memcpy(mxGetPr(x_L), x_Ldata, 2*1*sizeof(double));
  memcpy(mxGetPr(x_U), x_Udata, 2*1*sizeof(double));

  /* Call the Solveqp-routine */
  x = f = NULL;

  printf("Calling the solve routine...\n");

#ifdef MCC4
  mlfSolveqp(2, &x, &f, F, c, A, b_L, b_U, x_L, x_U);
#else
  x = mlfSolveqp(&f, F, c, A, b_L, b_U, x_L, x_U);
#endif

  if(x)
  {
    if(mxGetN(x)*mxGetM(x) >= 2)
    {
      printf("\nSoulution vector x:\n");
      printf("  %f  %f", mxGetPr(x)[0], mxGetPr(x)[1]);
    }
    mxDestroyArray(x);
  }
  if(f)
  {
    if(!mxIsEmpty(f))
    {
      printf("\n\nObjective function value f:\n");
      printf("  %f", mxGetPr(f)[0]);
    }
    mxDestroyArray(f);
  }
  
  printf("\n\nTerminating stand alone library...\n");
  libcplexqpTerminate();

  mxDestroyArray(F);
  mxDestroyArray(c);
  mxDestroyArray(A);
  mxDestroyArray(b_L);
  mxDestroyArray(b_U);
  mxDestroyArray(x_L);
  mxDestroyArray(x_U);

#ifdef MCC4
  printf("Terminating application...\n");
  mclTerminateApplication();
#endif

  return 0;
}

/*
 * MODIFICATION LOG:
 *
 * 041112 frhe  Written
 */
