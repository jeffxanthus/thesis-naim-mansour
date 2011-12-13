function o = angle(a)
% tomCmplx/angle - Argument of complex-valued object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

o = atan2(imag(a), real(a));

