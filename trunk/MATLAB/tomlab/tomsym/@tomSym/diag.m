function A = diag(B,d)
% tomSym/diag - Overloaded function
%
% A = DIAG(B,d) extracts diagonal of B specified by d.
% A = DIAG(V,d) creates a diagonal matrix with V on the d:th diagonal.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    d = 0;
end

if min(size(B))~=1
    % A = DIAG(B,d) extracts diagonal of B specified by d.
    if d==0 && size(B,1)==size(B,2)
        A = getdiag(B);
    else
        A = vec(lookup(B,diag(reshape(1:numel(B),size(B)),d)));
    end
else
    % A = DIAG(V,d) creates a diagonal matrix with V on the d:th diagonal.
    if d==0
        A = setdiag(B);
    else
        if d>0
            ii = 1:length(B);
            jj = (1:length(B))+d;
            n = length(B)+d;
        else
            ii = (1:length(B))-d;
            jj = 1:length(B);
            n = length(B)-d;
        end
        A = sparse(ii,jj,B,n,n);
    end
end
   
        
