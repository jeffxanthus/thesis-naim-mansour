function y = isreal(x) %#ok
% tomSym/isreal - Overloaded function
%
% y = isreal(x) where x is a tomSym object returns 'true', because a tomSym
% cannot be complex in the current implementation.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

y = true;
