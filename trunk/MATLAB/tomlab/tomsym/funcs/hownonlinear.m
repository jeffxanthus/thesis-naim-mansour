function d = hownonlinear(f) %# ok
% hownonlinear - Evaluate how nonlinear a tomSym expression is
%
% d = hownonlinear(f) returns:
%
% d = 0, if f is a constant
% 
% d = 1, if f is linear in all variables that it contains
%
% d = 2, if f is quadratic (not necessarily convex)
%
% d > 2 if none of the above is true.
%
% Calling this function is usually much faster than coputing the derivative
% and checking whether that is a constant.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Values other than d=0 are only returned by the overloaded function.
d = 0;


