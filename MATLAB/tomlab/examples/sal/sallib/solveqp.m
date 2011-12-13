% solveqp
%
% m-function to simplify the interface between C and MATLAB. It
% takes QP problem data as input, creates and solves a problem and
% returns the results of interest.
%
% Fredrik Hellman, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.5.0$
% Written Nov 11, 2004.   Last modified Nov 11, 2004.
%
function [x, f] = solveqp(F, c, A, b_L, b_U, x_L, x_U)

Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U);

Result = tomRun('cplex', Prob, [], 0);

x = Result.x_k;
f = Result.f_k;

% MODIFICATION LOG:
%
% 041112 frhe  Written
