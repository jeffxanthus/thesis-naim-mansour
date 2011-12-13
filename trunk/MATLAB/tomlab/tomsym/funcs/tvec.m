function T = tvec(m,n)
% tvec - The vectorized transpose matrix
%
% The vectorized transpose matrix, TVEC(m,n), is the m*n by m*n
% permutation matrix having the property that:
% kron(B,A) = tvec(size(B,1),size(A,1))*kron(A,B)*tvec(size(A,2),size(B,1))
%
% inv(T) = T'

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Definition:
%
% T = spalloc(m*n,m*n,m*n);
% for i=1:m
%     for j=1:n
% 	    T = T + kron(E(m,n,i,j),E(n,m,j,i));
% 	end
% end
%
% function e = E(m,n,i,j)
%   e = spalloc(m,n,1);
%   e(i,j) = 1;

if(nargin<2)
    n = m(2);
    m = m(1);
end

idx  = reshape(1:(m*n),n,m);
T = speye(m*n);
T = T(:,vec(idx'));
