function A = setSymmetric(v)
% tomSym/setSymmetric - Create a symmetric matrix
%
% A = setSymmetric(v) creates an n-by-n symmetric matrix A from the
% vector v of length n*(n+1)/2, by filling the lower triangluar part,
% column first, and mirroring to the upper triangular part.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~(size(v,1)==1 || size(v,2)==1)
    error('v must be a vector');
end

n = -0.5 + 0.5*sqrt(1+8*length(v));
if n~=round(n)
    error('Illegal length for v');
end

% Simplify in the scalar case
if n==1
    A = v;
else
    A = tomSym(mfilename,n,n,v);
end
