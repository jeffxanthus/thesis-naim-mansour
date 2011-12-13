function o = log(z)
% tomCmplx/log - Logarithm of complex-valued object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

o = tomCmplx(log(abs(z)), angle(z));
