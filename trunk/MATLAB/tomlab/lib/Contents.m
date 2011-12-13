% TOMLAB - A Toolbox for Applied Optimization
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% tomlab\lib  Library for TOMLAB
%
% Contents     This file
%
% AkaikeIC     Compute AIC; Akaike's information criterion (approximated)
% aicc         Compute AICc; The bias corrected Akaike's information criterion (approximated)
% AppRowQR     Append row to QR decomposition
% a2frstar     Convert node-arc A matrix to Forward-reverse star representation
% BoundInit    Initialize and check variables and bounds
% aVbuild      Compute alphaV=step length values for external solvers.
% backsubL     Solves lower triangular linear system; forward  substitution
% backsubR     Solves upper triangular linear system; backward substitution
% bic          Compute BIC; The Bayesian information criterion (approximated)
% binbin2lin   Converts binary products by adding linear constraints
% bincont2lin  Converts binary-integer/continuous products by adding linear constraints
% bmiran       Generator of random BMI (and quadratic MI) problems
% BoundInit    Initialization of bounds
% checkMAD     Check if MAD TB is installed correctly
% checkType    Check if input type (probType) same as (strType)
% checkx0      Check if x0 inside bounds, otherwise call pbuild.
% checkuP      Check if user parameters uP is properly set for current problem
% ComputeQR    Compute QR factorization
% con_fm       Compute merit function theta(x_k) used in constrained optim
% con_gm       Compute gradient of merit function theta(x_k) in con_fm
% consviolation  Compute constraint violations, return L1-norm
% convrate     Compute estimate of convergence rate
% conzero      Used to find x and y so constraint index equals zero.
% cpTransf     Transform general convex program to other forms
% CreateProbQP Create a structure for a subproblem of type QP,LP,FP or DLP.
% DefPar       Define structure parameter
% DeleteQR     Delete one column from the QR factorization
% eloLine      Line search algorithm for efficient local optimization
% endSolve     Bookkeeping, timings, compute search directions for optimization
% estBestHessian  Estimates the best hessian by selecting the step length
% estConsPattern  Estimates the ConsPattern
% estd2cPattern   Estimates the d2cPattern
% estHessPattern  Estimates the HessPattern
% estJacPattern   Estimates the Jacobian Pattern
% exp_eq       Computes a utility function needed by exp_q for zero finding
% exp_q        Find starting values for exponential parameters lambda
% expFitW      Retrive lambda and alpha from Prob.ExpFit.
% expGet1      Retrive five parameters from Prob.ExpFit.
% expGet2      Retrive four parameters from Prob.ExpFit.
% expGetLa     Retrive lambda and alpha from Prob.ExpFit.
% expInit      Find initial values for EXPFIT problem (Prob.ExpFit)
% ExpLinCon    Small code to generate linear constraints for ExpFit problem
% expLS        Solve NNLS problem for linear weights given exponential params
% expSet1      Set five parameters from Prob.ExpFit.
% expSet2      Set four parameters from Prob.ExpFit.
% expSetLa     Set lambda and alpha from Prob.ExpFit.
% fdng         Compute a Finite Numerical Difference Gradient (FNDG).
% fdng2        Compute a Finite Numerical Difference Gradient using splines
% fdng3        Compute a Fin. Numerical Diff. Gradient using complex vars
% fdnJ2        Compute a Finite Numerical Difference Jacobian using splines
% FDJac        Fin.Diff. Jacobian of residual or constraints using FD alg.
% FDJac2       Fin.Diff. Jacobian of residual or constraints using splines
% FDJac3       Fin.Diff. Jacobian of resid or constr. using complex vars
% FDcHess      Numerical approximation of the nonlinear constraints Hessian matrix.
% FDHess       Fin.Diff. Hessian using Gill et. al FD algorithm.
% FDHess2      Fin.Diff. Hessian using splines
% glbPrint     Print results during glbFast run
% glcPrint     Print results during glcFast run
% globalSave   Save global variables. Used in recursive optimization
% globalGet    Retrive global variables. Used in recursive optimization
% goal_c       Goal attainment constraints and residuals
% goal_dc      Goal attainment gradient of the constraints and jacobian of residuals
% goptions     Setup default values for optimization parameters and display
%              parameter values. Generalized version of OPTIM TB foptions.m
% HxFunc       Compute H*x for quadratic problems
%
% iniSolve     Initialization routine for optimization bookkeeping and timings
%              Init vars building solver search directions and step lengths
% inputSet     Prompt user in Window or GUI for integer from set of values
% inputR       Prompt user in Window or GUI for real input value in interval
% inputV       Prompt user for n-vector in interval. May duplicate 1 value
% InsertQR     Insert one column into the QR factorization
%
% ksDef        Index vectors and counters to be used in the Schittkowski routines
% ktr_d2L      KNITRO computes the Hessian
%
% L1_c         L1 constraints c and r
% L1_dc        L1 gradient of the constraints c and Jacobian
% L1_f         L1 objective function f
% L1_g         L1 gradient of the objective function f
% L1_H         L1 2nd derivative of the objective function f
% L2_c         L2 constraints c and r
% L2_dc        L2 gradient of the constraints c and Jacobian
% L2_f         L2 objective function f
% L2_g         L2 gradient of the objective function f
% L2_H         L2 2nd derivative of the objective function f
%
% LagMult      Compute Lagrange multiplier estimates
% LinearConstr Initialize and check linear constraints
% LineParamDef Define structure LineParam with line search parameters
% LineParamSet Set any parameters in structure LineParam that are undefined
% LinePlot     Plot along a line (the line search problem)
% LineSearch   Modified version of Fletcher's line search algorithm
% lls_r        Compute residual in Linear Least Squares Problem
% lls_J        Compute Jacobian in Linear Least Squares Problem
% lls_H        Compute zero Hessian in Linear Least Squares Problem
% lp_f         Define objective function for LP, c'*x (used by NLP solvers)
% lp_g         Define gradient vector for LP, c (used by NLP solvers)
% lp_H         Define Hessian matrix for LP, 0-matrix (used by NLP solvers)
% ls_f         Compute function value f(x) = 0.5 r(x)^T * r(x)
% ls_g         Compute gradient g(x) = J^T * r
% ls_H         Compute Hessian H(x) = J^T * J + d2r (2nd der part)
% ls_rJ        Computes both residual r(x) and J(x). Used with NLSSOL
% ls_rJS       Computes residual r(x) and sparse/dense J(x), used with OPTTB 2.x
% LSweights    Example routine, compute weights for NLLS, called by nlp_r
%
% mkbound      Make all bound variables defined and with same length
%
% mpec_f       Objective for MPEC problem
% mpec_g       Gradient for MPEC problem
% mpec_H       Hessian for MPEC problem
% mpec_c       Constraints for MPEC problem
% mpec_dc      Jacobian for MPEC problem
% mpec_d2c     Constraint Hessian for MPEC problem
%
% mPrint       Print matrix, format: NAME(i,:) a(i,1) a(i,2) ... a(i,n)
%
% nlp_c        Compute c using  Prob.FUNCS.c name of constraints
% nlp_cdc      Compute constraint c and derivatives dc using p_c and p_dc
% nlp_cdceq    Constraints c and derivatives dc for OPTIM TB 2.x
% nlp_cdcS     Compute c and dc (sparse) using p_c and p_dc (MINOS)
% nlp_cF       Compute  c with both linear and nonlinear constraints.
%              nlp_cF and nlp_dcF used by OPTIM TB-routines.
% nlp_cX       Compute  c with both linear and nonlinear constraints.
% nlp_d2c      Compute d2c using Prob.FUNCS.d2c name
% nlp_d2r      Compute d2r using Prob.FUNCS.d2r name
% nlp_dc       Compute dc using Prob.FUNCS.dc name of Jacobian or num.diff
% nlp_dcF      Compute dc with both linear and nonlinear constraints.
% nlp_dcX      Compute dc with both linear and nonlinear constraints.
% nlp_f        Compute f using  Prob.FUNCS.f name of function
% nlp_fc       Compute f and constraint c using  p_f and p_c (OPTIM)
% nlp_fg       Compute f and gradient g using  p_f and p_g.(MINOS)
% nlp_fgH      f,g and H for OPTTB 2.x
% nlp_g        Compute g using  Prob.FUNCS.g name of gradient or num.diff.
% nlp_gdc      Compute g and dc using p_g and p_dc (OPTIM)
% nlp_H        Compute H using  Prob.FUNCS.H name of Hessian  or num.diff.
% nlp_J        Compute J using  Prob.FUNCS.J name of Jacobian matrix func
% nlp_JT       Compute transposed Jacobian calling  p_J (OPTIM)
% nlp_r        Compute r using  Prob.FUNCS.r name of residual function
% nnls1        Solve sequence of nonnegative linear least squares problems
% NonlinConstr Initialize and check nonlinear constraints
%
% cpt_c        CONOPT to evaluate the nonlinear constraints
% cpt_d2L      CONOPT to computes the Hessian
% cpt_dc       CONOPT to computes Jacobian of the nonlinear constraints
%
% dffunc       Used as callback from Schittkowski routines
% dfgrad       Used as callback from Schittkowski routines
%
% opt15Run     Driver to call solvers from Optimization TB 1.5
% opt20Run     Driver to call solvers from Optimization TB 2.x
% optim_rJ     Compute residual r(x) and Jacobian in OPTTB 2.x interface
% optim_J      Return Jacobian J(x) in OPTTB 2.x interface
% optParamDef  Define structure optParam with optimization parameter fields
% optParamSet  Set any parameters in structure optParam that are undefined
%
% pbuild       Build solver search directions from iteration points
% plotInit     Set global structure plotData
% plotMenu     Menu for the plot options
% plotUtil     Plot utility routine for all plots in Tomlab
% preSolve     Presolve analysis on linear constraints, algorithm by Gondzio
% PrintMatrix  Print matrix with row and column labels
% PrintResult  Print optimization results
% PrintTSP     Print travelling saleman problem (TSP) results
% ProbCheck    Check that all values are set in structure Prob
% ProfDef      Initialization of input structure Prob.
% probInit     Interface routine for problem definition setup using ???_prob.
%
% qp_f         Compute QP function value
% qp_g         Compute QP gradient value
% qp_H         Compute QP Hessian value
% ResultDef    Define output result structure.
%
% Sep_rJ       Compute projections in separable NLLS algorithm
% strmenu      Converts string matrix to string. Run Matlab menu.
%
% smoothmax    Smoothing of max of variables as objective
% smoothmin    Smoothing of min of variables as objective
%
% WarmDefSOL   Move fields from Result.SOL to Prob.SOL for warm starts.
% vPrint       Print vector in rows, format: NAME(i1:in) v(i1) v(i2) ... v(in)
% xnargin      Utility to cope with a bug in Matlab 5.1. Calls nargin.
% z2frstar     Convert table of arcs (and costs) to Forward-reverse star
%
% expProbSet   Set exp values of input structure Prob
% expVarDef    Set initial values of exp variables
%
% GUI and menu utilities
% askAlgorithm.m
% askFlag         Common part of Menu code, reversing flag ask.
% askMethod       Ask for solver submethod choice
% askparam        Asks for one parameter in window or GUI. Used by xxx_prob.
% SolverAlgorithm
% SolverList      Gives list of solvers for a solvType
% SolverMethod    Return solver header and method menu for GUI and menus.
