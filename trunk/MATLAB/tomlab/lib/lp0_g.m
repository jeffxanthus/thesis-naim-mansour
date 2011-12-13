% function g = lp_g(x, Prob)
%
% Compute gradient to linear programming problem
%
% MINOS and SNOPT have explicit handling of linear term.
%
% Just put a 0 gradient here for the nonlinear part
%
% x      Point x where g(x) is evaluated 
% Prob   Problem structure
% g      Gradient, g(x).  g(x) = c;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 7, 2000.   Last modified Nov 7, 2000.

function g = lp0_g(x, Prob)

% Evaluate linear gradient
g = zeros(length(x),1);