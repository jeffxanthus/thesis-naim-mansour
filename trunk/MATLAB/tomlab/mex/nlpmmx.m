% NLPMMX Nonlinear Min-Max Matlab Solver
% -----------------------------------------------------------------------
%
%   NLPMMX solves the following nonlinear constrained min-max programming 
%   problem:
%
%   minimize    max ( f_1(x) ... f_L(x) )
%      x
%   subject to  x_L <=   x    <= x_U
%                      c_j(x)  =  0  ,  j=1,...,mEQ
%                      c_j(x) >=  0  ,  j=mEQ+1,...,m
%
%   where,
%
%   f(x) is an 1 x L vector of nonlinear functions.
%   c(x) is an m x 1 vector of nonlinear constraints.
%
% -----------------------------------------------------------------------
%
%   IF RUNNING TOMLAB:
%
%   Correct input is generated in nlpmmxTL.m from Prob input.
%   nlpmmx is called using the driver routine tomRun, or nlpmmxTL
%
%   IF NOT RUNNING TOMLAB:
%
%           If TOMLAB calls nlpmmx, then Prob is the standard
%           TOMLAB problem structure, otherwise the user should set:
%
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
%           Prob.ConsPattern = []; 
%
%           If the problem is a LP or QP problem (H defined), then the user 
%           does not have to specify anything more in the structure.
%
%           For a general nonlinear objective, or nonlinear constraints 
%           names of four user written routines must be given:
%
%           funr, actual name stored in Prob.FUNCS.r, with syntax
%                [r]     = funr(x, Prob)
%
%           funJ, actual name stored in Prob.FUNCS.J, wirh syntax
%                [J]     = funJ(x, Prob)
%
%           func, actual name stored in Prob.FUNCS.c, with syntax 
%                [c]     = func(x, Prob)
%
%           fundc, acutal name stored in Prob.FUNCS.dc, with syntax
%                [dc]    = func(x,Prob)
%
%           NLPMMX is calling the TOMLAB routines nlresid.m and nlJac.m
%           in the callback, and they call funr, func and funJ, fundc, respectively.
%
%           The user could also write his own versions of the routines
%           nlres.m and nlJac.m and put them before in the path.
%
% ------------------------------------------------------------------------
%
% function...
% [x, Inform, Iter, cLamda, Obj, f_eval, g_eval] = nlpmmx ( ...
%     L, m, mEQ, x_0, x_L, x_U, options, PriLev, PrintFile, Prob);
%
% ------------------------------------------------------------------------
% INPUT: (At least the six parameters must be given)
%
% L         The number of functions in the objective residual
% m         The number of constraints
% mEQ       The number of equality constraints
% x_0       Initial vector (DENSE).
% x_L       Lower bounds on x, n x 1 vector (DENSE).
% x_U       Upper bounds on x, n x 1 vector (DENSE).
% options   Optional options vector, with the following options:
% def  index name     description
% 1e-14 (1)  acc      Final termination accuracy.
% 1e-14 (2)  accqp    QP solver tolerance.
% 0.0   (3)  ressize  Guess for approximate size of objective.
% 20    (4)  maxfun   Upper bound on function calls during line search.
% 100   (5)  maxit    Maximum number of outer iterations. One iteration
%                     corresponds to either one evaluation of the gradients
%                     or one solution of the quadratic sub-problem.
% 0     (6)  maxnm    Stack size for merit function values during 
%                     non-monotonous line-search.
% 0.0   (7)  tolnm    Relative bound for merit function value increase.
%                     Only used if first step line search is unsuccesful.
%
% PriLev    Print level in the PrintFile. 
%        0   No output
%        1   Only a final convergence analysis is given.
%        2   One line of intermediate results is printed in each
%            iteration.
%        3   More detailed information is printed in each iteration
%            step, e.g. variable, constraint and multiplier values.
%        4   In addition to 'IPRINT=3', merit function and steplength
%            values are displayed during the line search.
%
% PrintFile The name of the file to be used for NLPMMX solver information 
%           output.
%
% Prob      Must be a structure. No check is made in the MEX interface.
%
% ------------------------------------------------------------------------
% OUTPUT:
% x         Solution vector with decision variable values (n x 1 vector).
% Inform    Result of NLPMMX run.
%           0 = Optimal solution found
% Iter      Number of iterations.
% cLamda    Multipliers (dual solution vector) (m+n+n x 1 vector).
% Obj       Objective function value at optimum.
% f_eval    The number of objective and constraint function evaluations.
% g_eval    The number of objective and constraint gradient evaluations.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 29, 2009.    Last modified Aug 19, 2009.
%# mex

function ...
    [x, Inform, Iter, cLamda, Obj, f_eval, g_eval] = ...
    nlpmmx( L, m, mEQ, x_0, x_L, x_U, ...
    options, PriLev, PrintFile, Prob);

help nlpmmx;