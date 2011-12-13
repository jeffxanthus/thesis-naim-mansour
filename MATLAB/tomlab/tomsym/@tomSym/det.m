function y = det(a)
% tomSym/det - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if size(a,1)~=size(a,2)
    error('Matrix must be square');
end

% Simplify if the argument is scalar or unity
if size(a,1)<=1
    y = a;
elseif isone(a)
    y = 1;
else
    y = tomSym(mfilename,1,1,a);
end
