% TOMLAB - Plug-In For Optimization TB
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% tomlab\optim
%
% ---------------------------------------------------------------------------
% Tomlab equivalents to Optimization Toolbox solvers 
% ---------------------------------------------------------------------------
%
% Contents     This file
%
% ---------------------------------------------------------------------------
% TOMLAB compatible interfaces are available for the following
% Optimization Toolbox routines:
% ---------------------------------------------------------------------------
%
% bintprog      Binary linear programming
% fgoalattain   Multidimensional multi-objective goal attainment
% fmincon       Constrained minimum of several variables
% fminimax      Minimax solution of a function of several variables
% fminunc       Multidimensional unconstrained nonlinear minimization
%               The Tomlab version handles simple bounds
% fseminf       Semiinfinite constrained optimization
%
% fsolve      - Nonlinear system of equations solve (function solve).
%               Use a standard nonlinear least squares formulation
%               in Tomlab to solve nonlinear system of equations
%
% linprog       Linear programming
%
% lsqcurvefit   Solves nonlinear least squares problems
% lsqlin        Linear least squares problem with linear constraints
% lsqnonlin     Solves nonlinear least squares problems
%
% quadprog      Quadratic programming
%
% ---------------------------------------------------------------------------
% Other optimization routines available in the general Matlab distribution
% ---------------------------------------------------------------------------
%
% fminbnd       Scalar bounded nonlinear function minimization
%               Tomlab Tfzero gives a faster, robust alternative
%
% fminsearch    Multidimensional unconstrained nonlinear minimization, 
%               by Nelder-Mead direct search method.
%
% fzero       - Scalar nonlinear zero finding.
%               Tomlab Tfzero gives a faster, robust alternative
%
% lsqnonneg   - Linear least squares with nonnegativity constraints.
%               Tomlab Tnnls gives a faster and more robust alternative
%               Tnnls also handles linear equalities
%
% lsqr          Large-scale Linear least squares and linear systems
%
% ---------------------------------------------------------------------------
% Tomlab has faster equivalents for some of the general optimization routines
% ---------------------------------------------------------------------------
%
% Tfzero        Tomlab Tfzero is a faster, robust alternative to fzero
%               and fminbnd
%
% Tlsqnoneng    An interface calling Tnnls with similar call as lsqnonneg
%
% Tlsqr         Solves large-scale Linear least squares and linear systems
%               with the original LSQR algorithm from Saunders and Paige.
%               Much faster and more robust, running with MEX file interface

%
% fminimax and fseminf (and fsolve) equivalents are not available.
% See Tomlab infSolve instead of fminimax
% See Tomlab NLLS routines instead of fsolve, to solve nonlinear equations
