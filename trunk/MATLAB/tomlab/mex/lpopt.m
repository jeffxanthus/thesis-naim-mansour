% LPOPT LP Matlab Solver
% -----------------------------------------------------------------------
%
%   lpopt solves the following linear programming problem (LP):
%
%        minimize  c'x    subject to:
%         x                            (  x )
%                               bl <=  ( Ax )  <= bu
%   where
%
%   A is an m x n dense matrix (linear constraints).
%   c is an n x 1 dense vector of linear objective coefficients.
%
%   bl, bu have dimension n+m
%
%   If isempty(c) (or all(c==0), then a feasible point problem is solved (FP).
%   Otherwise a standard linear programming problem is solved (LP).
%
% -----------------------------------------------------------------------
%
%   IF RUNNING TOMLAB:
%
%   bl, bu and other input are generated in lpoptTL.m from Prob input.
%   lpopt is called using the the driver routine tomRun, or calling lpoptTL
%
% ------------------------------------------------------------------------
%
% function [Inform, Iter, iState, Ax, cLamda, Obj, x] = lpopt( ...
%           A, bl, bu, c, Warm, x, iState, ...
%           SpecsFile, PrintFile, SummFile, PriLev, optPar );
%
% INPUT: (At least the first 3 parameters must be given)
%
% A         Constraint matrix, m x n (DENSE).
% bl        Lower bounds on (x,Ax), m+n x 1 vector (DENSE).
% bu        Upper bounds on (x,Ax), m+n x 1 vector (DENSE).
% c         Linear objective function cost coeffs, n x 1 (DENSE).
%           If length(c) < n, setting c(1:n)=0;
% Warm      If Warm > 0, then warm start, otherwise cold Start. Default 0.
%           If warm start, then x and iState must be set properly.
%           Normally the values from last call to lpopt are used.
% x         Initial estimate of solution vector x. (DENSE) If length(x) < n,
%           the rest of the elements in x are set to 0.
% iState    Working set (if Warm start) (n+m) x 1 (DENSE).
%           If length(iState) < n+m, setting iState(1:n+m)=0;
%
% iState(i)=0: Corresponding constraint not in the initial QP working set.
% iState(i)=1: Inequality constraint at its lower bound in QP working set.
% iState(i)=2: Inequality constraint at its upper bound in QP working set.
% iState(i)=3: Equality constraint in the initial QP working set,bl(i)==bu(i).
%
% SpecsFile Name of the OPTIONS File, see the TOMLAB /MINOS guide
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number 
%           of characters = 500
% PriLev    Print level in the lpopt MEX-interface.
%           = 0  Silent
%           = 1  Summary information
%           = 2  More detailed information
%           if isempty(PriLev), set as 0.
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           If length(optPar) < 62, lpopt sets the rest of the values to
%           missing value (-999).
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by LPOPT, unless the value is
% altered in the SPECS file. Refer to lpoptTL.m for information about optPar
% settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT:
% Inform    Result of LPOPT run.
%           0 = Optimal solution found
% Iter      Number of iterations.
% iState    Status of working set, se input description of iState.
% Ax        A*x.
% cLamda    Lagrangian multipliers (dual solution vector) (m x 1 vector).
% Obj       Objective function value at optimum.
% x         Solution vector with decision variable values (n x 1 vector).

% Kenneth Holmstrom, Tomlab Optimization AB, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization AB, Sweden. $Release: 4.8.0$
% Written July 12, 2000.     Last modified Jun 14, 2005.
%# mex

function [Inform, Iter, iState, Ax, cLamda, Obj, x] = lpopt( ...
          A, bl, bu, c, Warm, x, iState, ...
          SpecsFile, PrintFile, SummFile, PriLev, optPar )

help lpopt;