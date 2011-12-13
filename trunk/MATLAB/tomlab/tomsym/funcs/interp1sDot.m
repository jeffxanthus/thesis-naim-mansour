function yi = interp1sDot(x,y,xi,ndots)
% tomSym/interp1sDot - Derivative of interp1s
%
% yi = interp1sDot(x,y,xi,n) n:th derivative w.r.t xi of interp1s.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if size(x,2)~=1
    error('X must be a column vector');
end

if size(x,1)~=size(y,1)
    error('Input arguments x and y must have the same length.');
end

if nargin < 4
    ndots = 1;
end

if ndots==0
    yi = interp1s(x,y,xi);
else
    M = splineInterpolationMatrix(x,xi(:),ndots);
    if size(y,2)==1
        yi = reshape(M*y,size(xi,1),size(xi,2));
    else
        yi = M*y;
    end
end

