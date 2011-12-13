function y = negativeSemidefinite(A)
% negativeSemidefinite - Test if a matrix is negative definite
%
% y = negativeSemidefinite(A) is equivaltent to y = positiveSemidefinite(-A)
%
% See: positiveSemidefinite

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = positiveSemidefinite(-A);
