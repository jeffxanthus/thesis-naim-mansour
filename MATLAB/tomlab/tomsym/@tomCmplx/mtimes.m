function o = mtimes(a,b)
% tomCmplx/mtimes - Matrix multiplication

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-18 by rutquist for TOMLAB release 7.7

if size(a,1)==1 && size(b,2)==1 && ( isequal(a',b) || isequal(a,b') );
    o = real(a)*real(b)  - imag(a)*imag(b);
else
    o = tomCmplx(real(a)*real(b) - imag(a)*imag(b), ...
        real(a)*imag(b) + imag(a)*real(b));
end

