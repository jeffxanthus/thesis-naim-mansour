% LSSOL QP/LP/LS Matlab Solver
% -----------------------------------------------------------------------
%
% lssol solves the following quadratic programming problem (QP):
%
%      minimize  0.5*x'*H*x + c'x    subject to:
%         x                                   (  x )
%                                      bl <=  ( Ax )  <= bu
%
% or the least squares problem
%
%      minimize  0.5*(y-H*x)' * (y-H*x) + c'x    subject to:
%         x                                   (  x )
%                                      bl <=  ( Ax )  <= bu
% where,
%
% H is an n x n or m x n dense matrix (or empty, then a LP is solved).
% A is an nclin x n dense matrix (linear constraints) (or empty).
% c is an n x 1 vector of linear objective coefficients (or empty).
% y is an m x 1 vector of data observations to be fitted (or empty).
%
% bl, bu have dimension n+nclin.
% Upper is > 0 ==> H is upper trapezodial.
%
% Dependent on the input of H, y, c and Upper, any of the following ten
% problems are solved. If m == length(y) > 0, then a LS problem is assumed.
%
%     H      y      c      Upper
%   empty  empty   empty    --    Feasible point problem        (FP)
%   empty  empty   n x 1    --    Linear programming problem    (LP)
%   n x n  empty   empty     0    Quadratic programming problem (QP1)
%   n x n  empty   n x 1     0    Quadratic programming problem (QP2)
%   n x n  empty   empty     1    Quadratic programming problem (QP3)
%   n x n  empty   n x 1     1    Quadratic programming problem (QP4)
%   m x n    m     empty     0    Least squares problem         (LS1)
%   m x n    m     n x 1     0    Least squares problem         (LS2)
%   m x n    m     empty     1    Least squares problem         (LS3)
%   m x n    m     n x 1     1    Least squares problem         (LS4)
%
% -----------------------------------------------------------------------
%
% IF RUNNING TOMLAB:
%
% bl and bu are generated in lssolTL.m from Prob input.
% lssol is called using the driver routine tomRun or calling lssolTL.
%
% ------------------------------------------------------------------------
%
% function [x, Inform, iState, cLamda, Iter, fObj, r, kx] = ...
%     lssol(A, bl, bu, c, x, optPar, H, y, Warm, iState, Upper, kx, ...
%            SpecsFile, PrintFile, SummFile, PriLev, ProbName );
%
% INPUT: (At least the first 4 parameters must be given)
%
% A         Constraint matrix, nclin x n (DENSE).
% bl        Lower bounds on (x,Ax), nclin+n x 1 vector (DENSE).
% bu        Upper bounds on (x,Ax), nclin+n x 1 vector (DENSE).
% c         Linear objective function cost coeffs, n x 1 (DENSE).
%           If isempty(c), setting c(1:n)=0;
% x         Initial estimate of solution vector x. If isempty(x), x(1:n)=0.
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           If length(optPar) < 62, lssol sets the rest of the values to
%           missing value (-999).
% H         Matrix H in quadratic part of objective function (DENSE).
%           H is either empty (FP,LP), or quadratic n x n (QP), or m x n (LS).
% y         Data vector of length m for LS problem, otherwise empty.
% Warm      If Warm > 0, then warm start, otherwise cold Start. Default 0.
%           If warm start, then x and iState must be set properly.
%           Normally the values from last call to lssol are used.
% iState    Working set (if Warm start) (nclin+m) x 1 (DENSE).
%           If length(iState) < nclin+m, setting iState(1:n+m)=0 & Warm=0
%
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
%
% Upper     If > 0, then H must be an upper-trapezoidal matrix.
% kx        Order of the n columns in A for the QP3,QP4,LS3 or LS4 problem.
% SpecsFile Name of the OPTIONS File, see TOMLAB /SOL guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number 
%           of characters = 500
% PriLev    Print level in the lssol MEX-interface.
%           = 0  Silent
%           = 10 Dimensions are printed
%           if isempty(PriLev), set as 0.
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by LSSOL, unless the value
% is altered in the SPECS file. Refer to lssolTL.m for information about 
% optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT:
% x         Solution vector with decision variable values (n x 1 vector).
% Inform    Result of LSSOL run.
%           0 = x is a strong local optimum
%           1 = x is a weak local optimum (nonunique)
%           2 = the solutions seems to be unbounded
%           3 = No feasible point found
%           4 = Maximal number of iterations reached
%           5 = 50 changes of working set without change in x, cycling?
%           6 = An input parameter is invalid
% iState    Status of working set, se input description of iState.
% cLamda    Lagrangian multipliers (dual solution vector) (m x 1 vector).
% Iter      Number of iterations.
% Ax        A*x.
% fObj      Value of objective function if feasible, or sum of infeasibilities.
% r         Residual for LS problem.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.4.0$
% Written Nov 2, 2000.    Last modified Jul 14, 2006.
%# mex

function [x, Inform, iState, cLamda, Iter, fObj, r, kx] = ...
    lssol( A, bl, bu, c, x, optPar, H, y, Warm, iState, Upper, kx, ...
           SpecsFile, PrintFile, SummFile, PriLev, ProbName )

help lssol;