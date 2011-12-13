% function d2c = maked2cPatt(C)
%
% maked2cPatt computes the d2c matrix pattern using the constraint Jacobian
% pattern Prob.ConsPattern (in input C)
%
% INPUT:
% C    A 0-1 m by n-matrix, sparse or dense, with the constraint pattern
%      see the description of Prob.ConsPattern
%
% OUTPUT:
% d2c  A 0-1 n by n-matrix, sparse, with the sparsity pattern of the symmetric
%      second derivative matrix of the weighted sum of the constraints

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 11, 2004.  Last modified Apr 11, 2004.

function d2c = maked2cPatt(C)

[m,n]   = size(C);
[ix,iy] = find(C);

d2c     = sparse(n,n);

for i=1:n
    iz = find(i==iy);
    if ~isempty(iz)
       ixz = ix(iz);
       if length(ixz) == 1
          ic = find(C(ixz,:));
       else
          ic = find(any(C(ixz,:)));
       end
       d2c(ic,i) = 1;
    end
end

% MODIFICATION LOG
%
% 040411  hkh  Written