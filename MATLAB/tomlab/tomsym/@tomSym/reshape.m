function y = reshape(x,m,n)
% tomSym/reshape - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-03 by rutquist for TOMLAB release 7.7

if nargin<3
    if numel(m)==2
        n = m(2);
        m = m(1);
    elseif numel(m)<2
        error('Size vector must have at least two elements.');
    else
        warning('tomSym:reshape:ND',...
            'tomSym does not yet support N-D matrices. Result has been flattened to 2D.');
        n = prod(m(2:end));
        m = m(1);
    end
end

if numel(x)~=m*n
    error('The number of elements must not change.');
end

if size(x,1)==m && size(x,2)==n
    y = x;
elseif tomCmp(x,'vec') || tomCmp(x,'reshape') 
    y = reshape(operand(1,x),m,n);
elseif tomCmp(x,'ctranspose') && any(size(x)==1)
    y = reshape(operand(1,x),m,n);
elseif tomCmp(x,'sparse');
    [ii,jj] = ind2sub([m,n],sub2ind(size(x),operand(1,x),operand(2,x)));
    y = sparse(ii,jj,operand(3,x),m,n);
elseif n==1
    y = vec(x);
elseif m==1 && size(x,2)==1
    y = x';
else
    y = tomSym(mfilename,m,n,x,m,n);
end
