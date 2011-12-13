% TOMLAB - Costly Global Optimization
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% -------------------------------------------------------------------
% The CGO solvers, in directory cgo, mainly written in Matlab code.
% -------------------------------------------------------------------
% rbfSolve     Costly Global Optimization (in /CGO), RBF interpolation
% ego          Efficient Global Optimization (EGO) ( in /CGO)
%
% CGOrun       Runs the predefined global optimization problems in e.g. glb_prob.m
% cgoplot      Utility to make plots in TOMLAB /CGO
%
% gn_f         Compute gn_f using tomsol dll
% rbf_c        does rescaling of x before calling the constraint routine
% rbf_d2c      does rescaling of x before calling the second order constraint
% rbf_dc       does rescaling of x before calling the constraint Jacobian routine
% sn_f         Compute s_n(y) using tomsol dll
% sn_g         Computes gradient of s_n(y)
%
% dace_f       Computes the Likelihood function for the EGO algorithm.
% daceInit     Finds "space-filling" grid (EGO).
% ego_c        Constraint on the condition number of the correlation matrix
% ego_cc       Does rescaling for constraints
% ego_f        Computes the expected improvement for the EGO algorithm
