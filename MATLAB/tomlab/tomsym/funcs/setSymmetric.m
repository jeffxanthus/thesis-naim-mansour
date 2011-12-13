function A = setSymmetric(v)
% setSymmetric - Create a symmetric matrix
%
% A = setSymmetric(v) creates an n-by-n symmetric matrix A from the
% vector v of length n*(n+1)/2, by filling the lower triangluar part,
% column first, and mirroring to the upper triangular part.
%
% A = setSymmetric(A) where A is a symmetic matrix, returns A unchanged.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-03-04 by rutquist for TOMLAB release 7.7

if size(v,1) == size(v,2) && isequalwithequalnans(v,v')
    A = v;
    return
end

v = v(:);

n = -0.5 + 0.5*sqrt(1+8*length(v));
if n~=round(n)
    error('Illegal length for v');
end

if issparse(v)
    A = spalloc(n,n,2*nnz(v));
elseif islogical(v)
    A = false(n,n);
else
    A = zeros(n,n);
end

j = 1;
for i=1:n
    A(i,i:n)   = v(j:j+n-i)';
    A(i+1:n,i) = v(j+1:j+n-i);
    j = j+n-i+1;
end
