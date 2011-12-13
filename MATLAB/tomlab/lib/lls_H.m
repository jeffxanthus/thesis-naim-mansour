% function H = lls_H(x, Prob)
%
% Compute zero Hessian in Linear Least Squares Problem
%
% x      Parameter vector x 
% Prob   Problem structure
%        ||Cx-d|| is minimized where
%        C = Prob.LS.C
%        d = Prob.LS.y
%        The Jacobian is the constant matrix C
%        The Hessian is the zero matrix

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 5, 2000.      Last modified Nov 5, 2000.

function H = lls_H(x, Prob)

n = length(x);
H = sparse(n,n);