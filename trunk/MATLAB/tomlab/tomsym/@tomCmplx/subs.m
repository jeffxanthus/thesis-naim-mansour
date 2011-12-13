function o = subs(z,varargin)
% tomCmplx/subs - Substitute symbols

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

x = subs(real(z),varargin{:});
y = subs(imag(z),varargin{:});

if isnumeric(x) && isnumeric(y)
    o = complex(x,y);
else
    o = tomCmplx(x,y);
end
