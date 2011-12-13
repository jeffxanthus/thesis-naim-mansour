function y = DiracDelta(x,order)
% tomSym/DiracDelta - Placeholder derivative of a discontinuous function.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    order = 1;
end

y = tomSym(mfilename,size(x,1),size(x,2),x,order);
