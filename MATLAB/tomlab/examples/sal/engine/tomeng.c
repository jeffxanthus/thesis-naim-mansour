/*
 *	tomeng.c
 *
 *      Exemplifies how to call Matlab and Tomlab from a stand alone
 *      application using the Matlab engine.
 *
 *      Copyright (c) 2005 by Tomlab Optimization, Inc, e-mail: tomlab@tomlab.biz
 *      Written: Apr 5, 2005.   Last modified: Apr 5, 2005.
 *      
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"
#define  BUFSIZE 256

int main(int argv, char *argc[])
{
  Engine *ep;
  mxArray *Fmat, *cmat, *Amat, *b_Lmat, *b_Umat, 
    *x_Lmat, *x_Umat, *f_kmat, *x_kmat, *Resultmat, *Probmat;

  char *tomlabdir;
  char buffer[BUFSIZE];
  char *solver;

  /* Problem data */
  double F[9]   = {3, -1, 0, -1, 2, -1, 0, -1, 1};
  double c[3]   = {1, 1, 1};
  double A[3]   = {1, 2, 1};
  double b_L[1] = {4};
  double b_U[1] = {4};
  double x_U[3] = { 1e20,  1e20,  1e20};
  double x_L[3] = {-1e20, -1e20, -1e20};

  int status;

  /*
   * This program is called with one argument: the path to the tomlab directory
   */

  Fmat = cmat = Amat = b_Lmat = b_Umat = x_Lmat = x_Umat = f_kmat = x_kmat = Resultmat = Probmat = NULL;
  ep = NULL;
  status = EXIT_FAILURE;

  if(argv < 3)
  {
#ifdef PC
    printf("\nUsage: %s tomlab_directory tomlab_solver\nExample: %s c:\\tomlab cplex\n",
	   argc[0], argc[0]);
#else
    printf("\nUsage: %s tomlab_directory tomlab_solver\nExample: %s /usr/local/tomlab cplex\n",
	   argc[0], argc[0]);
#endif
    goto TERMINATE;
  }

  tomlabdir = argc[1];
  solver    = argc[2];

  /*
   * Start the MATLAB engine by calling engOpen.
   */
  printf("Open MATLAB engine... ");
  if (!(ep = engOpen("matlab -nojvm"))) {
    printf("Failure.\n");
    fprintf(stderr, "\nCan't start MATLAB engine.\n");
    goto TERMINATE;
  }
  printf("Done.\n");

  /*
   * First we want cd to the tomlab directory and call the
   * startup script.
   */
	
  {
    char *cdstr = "cd ";
    if(strlen(tomlabdir) >= BUFSIZE-strlen(cdstr))
    {
      fprintf(stderr, "\nPath string too long.\n");
      goto TERMINATE;
    }
    else
    {
      sprintf(buffer, "%s%s", cdstr, tomlabdir);
    }
  }

  /* Change directory */
  printf("Change work directory of MATLAB engine... ");
  engEvalString(ep, buffer);
  printf("Done.\n");

  /* Call startup */
  printf("Start up Tomlab... ");
  engEvalString(ep, "startup;");
  printf("Done.\n");
  
  /* We want to solve this problem: 
   *
   * min    2*x1^2 + 2*x2^2 + x3^2 - 2*x1*x2 - 2*x2*x3 + x1 + x2 + x3
   *
   * s.t.   x1 + 2*x2 + x3 = 4
   *        -inf < x < inf
   *
   * This could be setup as a Tomlab QP problem:
   *
   * F   = [ 3    -1     0
   *        -1     2    -1
   *         0    -1     1 ];
   * c   = [ 1     1     1 ]';
   * A   = [ 1     2     1 ];
   * b_L = [ 4 ];
   * b_U = [ 4 ];
   * x_L = [ -inf -inf -inf]';
   * x_U = [  inf  inf  inf]';
   */

  /* Create the matrices dexcribing the problem */
  printf("Create problem matrices... ");
  Fmat = mxCreateDoubleMatrix(3, 3, mxREAL);
  memcpy(mxGetPr(Fmat), F, 9*sizeof(double));
	
  cmat = mxCreateDoubleMatrix(3, 1, mxREAL);
  memcpy(mxGetPr(cmat), c, 3*sizeof(double));

  Amat = mxCreateDoubleMatrix(1, 3, mxREAL);
  memcpy(mxGetPr(Amat), A, 3*sizeof(double));

  b_Lmat = mxCreateDoubleMatrix(1, 1, mxREAL);
  memcpy(mxGetPr(b_Lmat), b_L, 1*sizeof(double));

  b_Umat = mxCreateDoubleMatrix(1, 1, mxREAL);
  memcpy(mxGetPr(b_Umat), b_U, 1*sizeof(double));

  x_Lmat = mxCreateDoubleMatrix(3, 1, mxREAL);
  memcpy(mxGetPr(x_Lmat), x_L, 3*sizeof(double));

  x_Umat = mxCreateDoubleMatrix(3, 1, mxREAL);
  memcpy(mxGetPr(x_Umat), x_U, 3*sizeof(double));

  /* Put the variables to the engine workspace */
  engPutVariable(ep, "F", Fmat);
  engPutVariable(ep, "c", cmat);
  engPutVariable(ep, "A", Amat);
  engPutVariable(ep, "b_L", b_Lmat);
  engPutVariable(ep, "b_U", b_Umat);
  engPutVariable(ep, "x_L", x_Lmat);
  engPutVariable(ep, "x_U", x_Umat);
  printf("Done.\n");

  /* Call qpAssign to create a Tomlab problem */
  printf("Assign problem... ");
  engEvalString(ep, "Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U);");
  if(Probmat = engGetVariable(ep, "Prob"))
  {
    printf("Done.\n");
    mxDestroyArray(Probmat);
  }
  else
  {
    printf("Failure.\n");
    fprintf(stderr, "\nCouldn't assign Tomlab problem. Check your Tomlab installation.\n");
    goto TERMINATE;
  }

  /* Have the solver to solve the problem */
  printf("Solve problem... ");
  {
    char *runstr = "Result = tomRun('%s', Prob, 1);";
    if(BUFSIZE - strlen(runstr) - strlen(solver) < 0)
    {
      fprintf(stderr, "\nSolver string too long.\n");
      goto TERMINATE;
    }
    sprintf(buffer, runstr, solver);
  }
  engEvalString(ep, buffer);
  Resultmat = engGetVariable(ep, "Result");
  if(!Resultmat)
  {
    printf("Failure.\n");
    fprintf(stderr, "\nCouldn't solve problem. Check your Tomlab installation and Tomlab license.\n");
    goto TERMINATE;
  }
  printf("Done.\n");

  /* Get the function value and the optimal x vector */
  printf("Get results... ");
  f_kmat = mxGetField(Resultmat, 0, "f_k");
  if(!f_kmat)
  {
    printf("Failure.\n");
    fprintf(stderr, "\nCouldn't get the f_k value from the Result structure.\n");
    goto TERMINATE;
  }

  x_kmat = mxGetField(Resultmat, 0, "x_k");
  if(!x_kmat)
  {
    printf("Failure.\n");
    fprintf(stderr, "\nCouldn't get the x_k value from the Result structure.\n");
    goto TERMINATE;
  }
  printf("Done.\n");

  /* Print results */
  printf("\nRESULTS\n"
	 "--------------------------------------------------------------------\n");
  printf("Function value at optimal point: %f\n", *mxGetPr(f_kmat));
  printf("Optimal point:                   [%f  %f  %f]'\n", 
	 mxGetPr(x_kmat)[0], mxGetPr(x_kmat)[1], mxGetPr(x_kmat)[2]);
  printf("--------------------------------------------------------------------\n");

  status = EXIT_SUCCESS;
  printf("\nDone!\n");

  /* Deallocate Matlab matrices, if allocated */
 TERMINATE:
  if(Resultmat) mxDestroyArray(Resultmat);
  if(Fmat)      mxDestroyArray(Fmat);
  if(cmat)      mxDestroyArray(cmat);
  if(Amat)      mxDestroyArray(Amat);
  if(b_Lmat)    mxDestroyArray(b_Lmat);
  if(b_Umat)    mxDestroyArray(b_Umat);
  if(x_Lmat)    mxDestroyArray(x_Lmat);
  if(x_Umat)    mxDestroyArray(x_Umat);

  /* Close engine environment */
  if(ep)        engClose(ep);

  return status;
}







