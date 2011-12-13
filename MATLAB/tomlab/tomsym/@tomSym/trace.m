function y = trace(a)
% tomSym/trace - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if size(a,1)~=size(a,2)
    error('Matrix must be square');
end

% Simplify if the argument is unity or an inverse
y = tomSym(mfilename,1,1,a);
