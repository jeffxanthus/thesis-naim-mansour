function M = sparse(i,j,s,m,n)
% tomSym/sparse - Overloaded function
%
% M = SPARSE(i,j,s,m,n) generatres an m-by-n matrix M where
% M(i(k),j(k))=s(k), and the remaining elements of M are zero.

% SPARSE does not convert matrices internal to an object to sparse, but all
% constant matrices that tomSym creates are sparse from the start anyway.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-14 by rutquist for TOMLAB release 7.7

if nargin==1
    % Converting a tomSym to sparse has no effect.
    M = i;
    return;
end

if numel(i)~=numel(j) || numel(j)~=numel(s)
    error('I, J, and S must all have the same number of elements.');
end

if nargin<5 && isnumeric(j)
    n = max(j(:));
end

if nargin<4 && isnumeric(i)
    m = max(i(:));
end

if numel(s)>m*n
    error('Too many elements');
end

if (isnumeric(i) && (any(i(:)<1) || any(i(:)>m))) || ...
        (isnumeric(j) && (any(j(:)<1) || any(j(:)>n)))
    error('Index exceeds matrix dimensions.');
end

if tomCmp(s,'sparse')
    idx = sub2ind(size(s),operand(1,s),operand(2,s));
    i = i(idx);
    j = j(idx);
    s = operand(3,s);
end

if isempty(s)
    M = spzeros(m,n);
elseif m==1 && n==1
    M = s;
elseif m==n && isnumeric(i) && isnumeric(j) && numel(i)==m && ...
        numel(j)==n && all(i(:)==(1:m)') && all(j(:)==(1:n)')
    M = setdiag(s);
elseif numel(i)==m*n && isnumeric(i) && isnumeric(j) && ...
        all(sub2ind([m n], vec(i), vec(j))==vec(1:m*n))
    M = reshape(s,m,n);
else
    M = tomSym(mfilename,m,n,vec(i),vec(j),s,m,n);
end
