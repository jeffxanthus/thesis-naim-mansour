function yi = interp1sDot(x,y,xi,ndots)
% tomSym/interp1sDot - Derivative of interp1s, overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin < 4
    ndots = 1;
end

if size(x,2)~=1
    error('X must be a column vector');
end

if size(x,1)~=size(y,1)
    error('Input arguments x and y must have the same length.');
end

if ndots==0
    yi = interp1s(x,y,xi);
elseif isnumeric(x) && isnumeric(xi)
    M  = splineInterpolationMatrix(x,vec(xi),ndots);
    if size(y,2)==1
        yi = reshape(M*y,size(xi,1),size(xi,2));
    else
        yi = M*y;
    end
else
    if size(y,2)==1
        yi = tomSym(mfilename,size(xi,1),size(xi,2),x,y,xi,ndots);
    else
        yi = tomSym(mfilename,numel(xi),size(y,2),x,y,xi,ndots);
    end
end
