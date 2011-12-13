function yi = interp1l(x,y,xi)
% tomSym/interp1l - Linear interpolation
%
% yi = interp1p(x,y,xi) is equivalent to 
% yi = interp1(x,y,xi,'linear','extrap')
%
% This function is used internally by tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

yi = interp1(full(x),full(y),full(xi),'linear','extrap');
