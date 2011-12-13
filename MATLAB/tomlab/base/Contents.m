% TOMLAB - Base Module
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% ----------------------------------
% Driver routines
% ----------------------------------
% tomRun        General driver routine calling any Tomlab solver
% tomSolve      General subproblem driver routine calling any Tomlab solver
%               No check is made on input problem structure.
%               Global variables are saved and restored.
% tomRunMini    Silent driver routine for LP/QP subproblems
%
% -------------------------------------------
% Files in the BASE MODULE, in directory base
% -------------------------------------------
%
% tomlabInit    Tomlab global variable initialization (not necessary to run)
% tomlab.bib    Bibtex file with Tomlab references
% checkdll      Checks the licenses
%
% -----------------------------------------------
% The TOM solvers, in directory base, mainly written in Matlab code.
% -----------------------------------------------
% clsSolve      Constrained nonlinear least squares solver:
%               Active set strategy. As search step algorithm using:
%               Gauss-Newton with subspace minimization, Fletcher-Xu hybrid
%               method, Al-Baali-Fletcher hybrid method and Huschens method
% conSolve      SQP algorithm. Schittkowski with Augmented Lagrangian
%               Han-Powell (Quasi-Newton update)
% cutplane      Gomory's cutting plane algorithm for Mixed-Integer Programs
% DualSolve     Dual simplex algorithm with three selection rules.
% expSolve      Solve exponential fitting problems
% glbFast       glbSolve DIRECT algorithm in faster Fortran version
% glbSolve      Global Optimization algorithm DIRECT by Don Jones et.al.
% glcCluster    Hybrid of glcFast / cluster alg / local solver
% glcFast       glcSolve constrained DIRECT algorithm in faster Fortran version
% glcSolve      Constrained Mixed-Integer Global Optimization, Constr.DIRECT
% goalSolve     Multi-objective Goal Attainment
% infLinSolve   Linear minimax solver
% infSolve      Sparse constrained minimax solver
% L1LinSolve    Finds a linearly constrained L1 solution
% L1Solve       Sparse constrained L1 solver
% linRatSolve   Linear ratio solver
% lpSimplex     Simplex algorithm for general LP, structure input.
%               Solves Phase 1 and Phase 2 problems
% mipSolve      Branch & Bound for Mixed-Integer Programs (MIP)
% multiMin      Does local optimization from a set of starting points.
% nlpSolve      Filter SQP algorithm by Fletcher-Leyffer. Calls qpPhase1.
% qpSolve       Solve QP using active set method.
% slsSolve      Sparse Least Squares solver (constrained)
% sTrustr       A Structural Trust Region algorithm for unconstrained
%               optimization (Conn, Gould, Sartenaer, Toint). Calls itrr
% stepexp       Stepwise solve exponential fitting problems
% ucSolve       Unconstrained optimization solver, handling simple bounds
%               on the parameters. Algorithms: Newton, BFGS, Inverse BFGS
%
% -------------------------------------------------------------------
% Utilities to help defining problems in the Tomlab (TQ) format
% -------------------------------------------------------------------
% checkAssign   Check most of the inputs for assign routines below
% probAssign    Setup a Prob structure for a problem of certain problem type
% lpAssign      Define a Linear Programming problem
% lpconAssign   Define a Nonlinearly constrained LP problem
% qpAssign      Define a Quadratic Programming problem
% qpconAssign   Define a Nonlinearly constrained QP problem
% conAssign     Define a NonLinear Programming problem (constrained or not)
% mipAssign     Define a mixed-integer programming problem
% miqpAssign    Define a mixed-integer quadratic programming problem
% miqqAssign    Define a mixed-integer quadratic programming problem
%               with quadratic constraints
% clsAssign     Define a nonlinear least squares problem (constrained or not)
% llsAssign     Define a linear least squares problem (constrained or not)
% glcAssign     Define a global optimization problem (constrained or not)
% sdpAssign     Define a semidefinite program
% bmiAssign     Define a bilinear semidefinite program
% minlpAssign   Define a mixed-integer nonlinear (MINLP) program
% simAssign     Both the function and the constraints are computed
% expAssign     Define an exponential fitting problem
% tomFiles      Set names of the m-files to be used into structure Prob
%               Only needed if using the general probAssign routine
%
% BuildMPEC     Defines the structure for an MPEC problem
%
% keep_A        Keeps a set of linear constraints
% modify_b_L    Modifies the lower bounds on the linear constraints
% modify_b_U    Modifies the upper bounds on the linear constraints
% modify_c      Modifies the linear part for the objective function (LP, QP)
% modify_c_L    Modifies the lower bounds on the nonlinear constraints
% modify_c_U    Modifies the upper bounds on the nonlinear constraints
% modify_x_0    Modifies the starting point
% modify_x_L    Modifies the lower bounds on the decision variables
% modify_x_U    Modifies the upper bounds on the decision variables
% remove_A      Removes a set of linear constraints
% replace_A     Replaces the linear constraints
%
% CreateTomProb     Creates a file TomlabProblem.mat with the standard
%                   predefined test problems in Tomlab
%
% -----------------
% Testing utilities
% -----------------
% pretest       Test of presolve analysis, calls preSolve.m
% runtest       Run test of solver for sequence of problems
% systest       Run test of many solvers for sets of problems
% testtom       Script to run Tomlab system tests
