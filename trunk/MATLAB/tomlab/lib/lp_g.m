% function g = lp_g(x, Prob)
%
% Compute gradient to linear programming problem
%
% x      Point x where g(x) is evaluated
% Prob   Problem structure
% g      Gradient, g(x).  g(x) = c;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 5, 1998.     Last modified Sep 5, 1999.

function g = lp_g(x, Prob)
g = Prob.QP.c;