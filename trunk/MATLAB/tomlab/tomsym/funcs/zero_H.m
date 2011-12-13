function z = zero_H(x,Prob) %#ok
% zero_H - The Hessian matrix of zero_f
%
% z = zero_H(x,Prob) returns an all-zero numel(x)-by-numel(x) matrix.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

z = spalloc(numel(x),numel(x),0);
