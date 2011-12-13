function J = prodJ1(M,dim)
% tomSym/prodJ1 - derviative of prod
%
% J = prodJ1(M,dim) computes the Jacobian (derivative) matrix of the
% function call prod(M,dim)
%
% Each row of the returned matrix J corresponds to an element of
% prod(M,dim) and each column corresponds to an element of M.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

p = prod(M,dim);

sz = [numel(p), numel(M)];
szm = size(M);

M1 = M;
M1(M==0) = 1;

ix = reshape(1:sz(2),szm(1),szm(2));
if dim==1
    jx = repmat((1:szm(2)),szm(1),1);
    v  = repmat(p,szm(1),1)./M1;
    for i=find(p==0)
        j = find(M(:,i)==0);
        if length(j)>1
            v(j,i) = 0;
        else
            v(j,i) = prod(M1(:,i));
        end
    end
    J  = sparse(vec(jx),vec(ix),vec(v),sz(1),sz(2));
else
    jx = repmat((1:szm(1)),szm(2),1);
    v  = repmat(p,1,szm(2))./M1;
    for i=find(p'==0)
        j = find(M(i,:)==0);
        if length(j)>1
            v(i,j) = 0;
        else
            v(i,j) = prod(M1(i,:));
        end
    end
    J  = sparse(vec(jx'),vec(ix),vec(v),sz(1),sz(2));
end
