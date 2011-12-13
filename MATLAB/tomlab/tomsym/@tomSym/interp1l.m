function yi = interp1l(x,y,xi)
% tomSym/interp1l - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-26 by rutquist for TOMLAB release 7.7

if ~any(size(x)==1)
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

if length(x)==2
    % Simplify for the case with exactly two points.
    x1 = lookup(x,1);
    x2 = lookup(x,2);
    p = (vec(xi)-x1)./(x2-x1);
    yi = (1-p)*submatrix(y,1,':') + p*submatrix(y,2,':');
    if length(y)==numel(y)
        yi = reshape(yi,size(xi));
    end
    return;
end

if isequal(vec(x),vec(xi))
    if size(y,2)==1
        yi = reshape(y,size(xi));
    else
        yi = y;
    end
elseif isnumeric(x) && isnumeric(xi)
    M  = linearInterpolationMatrix(x,vec(xi));
    if size(y,2)==1
        yi = reshape(M*y,size(xi,1),size(xi,2));
    else
        yi = M*y;
    end
elseif isnumeric(x) && isnumeric(y)
    try
        pp = interp1(x,y,'linear','pp'); 
    catch
        pp = [];
    end
    if ~isstruct(pp);
        % TODO: Fix this to be compatible with Matlab 6.5
        warning('tomSym:interp1:OldMatlabVersion',...
            'This method of interpolation is not fully supported by tomSym under Matlab 6.5 or earlier');
        if size(y,2)==1
            yi = tomSym(mfilename,size(xi,1),size(xi,2),x,y,xi);
        else
            yi = tomSym(mfilename,numel(xi),size(y,2),x,y,xi);
        end
    else
        yi = ppval(pp,xi);
    end
else
    if size(y,2)==1
        yi = tomSym(mfilename,size(xi,1),size(xi,2),x,y,xi);
    else
        yi = tomSym(mfilename,numel(xi),size(y,2),x,y,xi);
    end
end
