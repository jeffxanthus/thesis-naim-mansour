% snoptA.m
%
% SNOPT 7.2-5 NLP Solver
%
% snoptA solves the following general NLP problem:
%
%   min Fobj(x)
%    x 
%
%   subject to     x_L <=   x  <= x_U
%                  b_L <= F(x) <= b_U
%
% where x in R^n is the vector of unknowns and F(x) in R^m is a vector of 
% nonlinear functions of x. 
%
% F(x) is given as F(x) = f(x) + A*x, i.e. a sum of a nonlinear component 
% f(x) in R^m  and a  linear component A*x, A in R^(m,n).
%
% The objective function Fobj(x) is one element of the F(x) vector. 
%
% -------------------------------------------------------------------------
% Call syntax: 
%
% [x_k, xstate, xmul, F, Fmul, Fstate, Inform, nS, nInf, sInf] = ...
%    snoptA( x_L, x_U, x_0, b_L, b_U, A, G, Func, objRow, optPar, ...
%            Warmstart, xState, fState, nS, Prob, ...
%            PrintFile, SummFile, PriLev );
%
% Inputs:
%
% x_L, x_U     Lower and upper variable bounds. For an unbounded variable x(i),
%              set x_L(i) = -Inf and/or x_U(i) = Inf.
%
% x_0          Starting point for snoptA. If given, x_0 must be dense and of 
%              the same length as x_L and x_U. 
%
% b_L, b_U     Lower and upper function bounds. For an unbounded function F(i),
%              set b_L(i) = -Inf and/or b_U(i) = Inf. In particular, the
%              objective row should generally be set as unlimited.
%
% A            The linear function matrix in F(x) = f(x) + A*x. Can be 
%              either sparse (recommended) or dense. 
%
% G            Nonzero element pattern of the nonlinear problem functions
%              f(x) - one row per function, one column per variable. Can be
%              sparse (recommended) or dense. 
%
% Func         Can be either a string containing the name of a valid M-file
%              (or MEX-file) that calculates the nonlinear function vector
%              and it's Jacobian or a handle (function_handle) to a similar
%              function. 
%
%              The syntax of the user function is described below.
%
% objRow       Index of the objective row. F(objRow) is taken as the
%              function Fobj(x) to minimize.
%
% optPar       Vector of SNOPT parameter values.
%
% Warmstart    Flag for Warm Start: default is 0, Cold Start. If set to 1,
%              xState, fState and nS must be given.
%
% xState       States of variables, used for Warm Start. 
%
% fState       States of functions, used for Warm Start.
%
% nS           Number of superbasic variables, used for Warm Start.
%
% Prob         A user settable field which is passed to the user function
%              callback.
%
% PrintFile    Name of SNOPT print file with detailed iteration
%              log. If no print file is desired, set to empty [] or ''
%
% SummFile     Name of SNOPT summary file with brief iteration
%              log. If no summary file is desired, set to empty [] or ''
%
% PriLev       Print Level in snoptA, controls the amount of information 
%              printed to the command window. The following values are 
%              recognized:
%
%              0 - Silent (default)
%              1 - Output similar to the Summary File
%              2 - Output similar to the Print File
%
% -------------------------------------------------------------------------
% Outputs:
%
% x_k          Final values of the variables
%              
% xstate       Final state of the variables as follows:
%              
%              xstate(j)   State of x(j)    Usual value of x(j)
%              
%                 0        nonbasic         x_L(j)
%                 1        nonbasic         x_U(j)
%                 2        superbasic       Between x_L(j) and x_U(j)
%                 3        basic            Between x_L(j) and x_U(j)
%              
% xmul         Vector of dual variables for the bound constraints
%              x_L <= x <= x_U
%              
% F            Final value of the vector of problem functions F(x) 
%              
% Fmul         Vector of dual variables (Lagrange multipliers) for
%              the general constraints b_L <= F(x) <= b_U
%              
% Fstate       Final state of the problem functions. 
%              
% Inform       Reports the result of the call to snoptA. 
%              
% nS           Final number of superbasic variables.
%              
% nInf         Final number of infeasibilities
%              
% sInf         Sum of the infeasibilities of constraints that lie
%              outside one of their bounds by more than optPar(???)
%              before the solution is unscaled. 
%
%
% User function syntax:
%
% The user function given in the Func input must adhere to the following
% syntax:
%
%
% function [F,G,mode] = UserFunc(x,status,Prob)
%
% INPUTS:
%
%  x         The decision variable vector
%
%  status    A 1x3 vector with the following information:
%
%  status(1) If snoptA calls the user function for the first time, this
%            element contains a 1 (one), zero (0) otherwise.
%
%  status(2) needF - flag indicating if F needs to be calculated.
%  status(3) needG - flag indicating if G needs to be calculated.
%
%  Prob      The Prob structure passed to snoptA. This can be used for
%            problem specific information. 
%
% OUTPUTS:
%
%  F         Vector of nonlinear function components. If status(2)==0,
%            an empty array [] may be returned in this output.  
%
%  G         Jacobian of nonlinear functions. One row per function, one
%            column per variable. If status(3)==0, an empty array [] may be
%            returned in this output. 
%
%  mode      Status flag from user routine. There are two values with special
%            significance to snoptA:
%
%            Returning mode=-1 tells snoptA that one or more problem functions
%            can not be evaluated at the current point x. snoptA linesearch will
%            shorten the step and try again. 
%
%            Returning any value mode < -1 instructs snoptA to stop the
%            optimization.
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written Feb 02, 2008.  Last modified Feb 2, 2008.
%# mex

help snoptA
