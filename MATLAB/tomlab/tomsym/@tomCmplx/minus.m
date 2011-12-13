function o = minus(a,b)
% tomCmplx/minus - Subtract object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-27 by rutquist for TOMLAB release 7.7

o = tomCmplx(real(a)-real(b), imag(a)-imag(b));

