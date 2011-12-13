function yi = interp1s(x,y,xi)
% tomSym/interp1s - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-22 by rutquist for TOMLAB release 7.7

if numel(x) ~= length(x)
    error('X must be a vector.');
end

if size(x,2)>1
    x = x';
end

if size(y,1)==1 && size(x,1)>1 && size(y,2)==size(x,1)
    y = y';
end

if size(x,1)~=size(y,1)
    error('Input arguments x and y must have the same length.');
end

if ( tomCmp(x,'smtimes') && (tomCmp(xi,'smtimes') || ...
        (tomCmp(xi,'mtimes') && numel(xi)==1)) || ...
        tomCmp(x,'plus') && tomCmp(xi,'plus') ) && ...
        isequal(operand(1,x),operand(1,xi))
    % Simplifications necessary for atPoints to work as expected.
    yi = interp1(operand(2,x),y,operand(2,xi),method);
    return;
end

if isnumeric(xi) && numel(xi)==1 && isnumeric(x)
    % Simplify for exact match.
    idx = find(x==xi);
    if ~isempty(idx)
        yi = lookup(y,idx);
        return
    end
end

if isequal(vec(x),vec(xi))
    if size(y,2)==1
        yi = reshape(y,size(xi));
    else
        yi = y;
    end
elseif isnumeric(x) && isnumeric(xi)
    M  = splineInterpolationMatrix(x,vec(xi));
    if size(y,2)==1
        yi = reshape(M*y,size(xi,1),size(xi,2));
    else
        yi = M*y;
    end
elseif isnumeric(x) && isnumeric(y)
    pp = spline(x',y');
    yi = ppval(pp,xi)';
else
    if size(y,2)==1
        yi = tomSym(mfilename,size(xi,1),size(xi,2),x,y,xi);
    else
        yi = tomSym(mfilename,numel(xi),size(y,2),x,y,xi);
    end
end
