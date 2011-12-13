% LSQL QP Matlab Solver
% -----------------------------------------------------------------------
%
%   LSQL solves the following convex quadratic programming problem (QP):
%
%   minimize    0.5 * x'* H * x + c' * x
%      x
%   subject to  A(j,1:n-L)*x1 - x2(j) +  b(j)   =  0  ,  j = 1,...,L
%               A(j,:) * x            +  b(j)   =  0  ,  j = L,...,mEQ
%               A(j,:) * x            +  b(j)  >=  0  ,  J = mEQ+1,...,m
%               x_L  <=  x1  <= x_U
%   where
%
%   x1 is the first n-L variables of x
%   x2 is the last L variables of x
%   H is an n x n dense positive definite matrix
%   A is an m x n dense matrix (linear constraints)
%     where the first L x n rows contain an L x L identity matrix
%   c is an n x 1 vector of linear objective coefficients.
%
%   If isempty(H) and isempty(c), then a feasible point problem is solved (FP)
%   If isempty(H), then a linear programming problem is solved (LP)
%   If isempty(c), then a quadratic programming problem is solved (QP1)
%   Otherwise a standard quadratic programming problem is solved (QP2)
%
% -----------------------------------------------------------------------
%
%   IF RUNNING TOMLAB:
%
%   Correct input is generated in lsqlTL.m from Prob input.
%   lsql is called using the driver routine tomRun, or lsqlTL
%
% ------------------------------------------------------------------------
%
% function [x, Inform, Iter, iState, Ax, cLamda, Obj] = lsql ( ...
%         H, A, b, mEQ, c, x_L, x_U);
%
% INPUT: (At least the first parameter H must be given)
%
% H         Matrix H in quadratic part of objective function (DENSE).
% A         Constraint matrix, m x n. (DENSE).
% b         Lower bounds on (Ax), m x 1 vector (DENSE).
% L         Number of special equality constraints
% mEQ       Number of equality constraint out of the m linear constraints.
% c         Linear objective function cost coeffs, n x 1 (DENSE).
%           If length(c) < n, setting cvec(1:n)=0;
% x_L       Lower bounds on x, n x 1 vector (DENSE). Default -1E20.
% x_U       Upper bounds on x, n x 1 vector (DENSE). Default  1E20.
%
% OUTPUT:
% x         Solution vector with decision variable values (n x 1 vector).
% Inform    Result of LSQL run.
%           0 = Optimal solution found
% Iter      Number of iterations.
% iState    Basis status of constraints + variables, (m + n x 1 vector)
%           State of variables: 0=nonbasic(on bl),1=nonbasic(on bu)
%                 2=superbasic (between bounds),3=basic (between bounds)
% Ax        A*x.
% cLamda    Lagrangian multipliers (dual solution vector) (m x 1 vector).
% Obj       Objective function value at optimum.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 1, 2000.    Last modified Feb 23, 2004.
%# mex

function [x, Inform, Iter, iState, Ax, cLamda, Obj] = lsql ( ...
        H, A, b, L, mEQ, c, x_L, x_U)

help lsql;