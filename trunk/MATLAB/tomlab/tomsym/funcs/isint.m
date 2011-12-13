function y = isint(x)
% isint - Determine if a number is integer
%
% y = isint(x) returns "true" if all elements of x equal round(x)
%
% isint should not be confused with Matlab's builtin "isinteger".
% isint(3) returns "true" while isinteger(3) returns false.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

y = ( all(x(:) == round(x(:))) );
