function o = times(a,b)
% tomCmplx/times - Elementwise multiplication

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

o = tomCmplx(real(a).*real(b) - imag(a).*imag(b), ...
    real(a).*imag(b) + imag(a).*real(b));
