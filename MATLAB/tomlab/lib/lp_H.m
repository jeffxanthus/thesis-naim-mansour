% function H = lp_H(x, Prob)
%
% Compute Hessian to linear problem, i.e. zero matrix
%
% x      Point x where H(x) is evaluated 
% Prob   Problem structure
% H      Hessian matrix, H(x) = 0, in f(x) = c'*x;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlan@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 5, 1998.    Last modified Jan 9, 2002.

function H = lp_H(x, Prob)
n=length(x);
H = sparse(n,n);