% Compute objective function value for linear programming problem
%
% MINOS and SNOPT have explicit handling of linear term.
%
% Just put 0 here for the nonlinear part
%
% function f = lp0_f(x, Prob)
%
% x      Point x where f(x) is evaluated 
% Prob   Problem structure
% f      Function value, f(x).  f(x) =  c'*x;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 7, 2000.   Last modified Nov 7, 2000.

function f = lp0_f(x, Prob)

% Evaluate linear objective function
f=0;