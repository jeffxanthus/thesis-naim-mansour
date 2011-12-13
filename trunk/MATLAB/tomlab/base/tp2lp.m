% tp2lp Convert Transportation LP to Linear Programming formulation
%       using supply, demand and distance matrix
%
% function [A,b,c] = tp2lp(s, d, C)
%
%                    m = M + N; n = M * N
%
%   Input:
%           s        Supply, N by 1
%           d        Demand, M by 1
%           C        Cost (Distance) matrix, M by N
%   Output:
%           A 	     Node-arc incidence matrix. A is m x n.
%           b 	     Supply/demand vector. Length m.
%           c 	     Cost vector. Length n.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Mar 1, 2003.   Last modified Mar 22, 2005.

function [A,b,c] = tp2lp(s, d, C)

m = length(s);
n = length(d);
mn = m*n;

i1 = ones(2*mn,1);
j1 = ones(2*mn,1);
v  = ones(2*mn,1);

k = 0;
l = 0;
for i = 1:m
    for j=1:n
        k = k + 1;
        l = l + 1;
        i1(l) = i;
        j1(l) = k;

        l = l + 1;
        i1(l) = m+j;
        j1(l) = k;
        v(l)  = -1;
    end
end
A = sparse(i1,j1,v);

b = [s(:);-d(:)];
C = C';
c = C(:);

% MODIFICATION LOG
%
% 030301 hkh Written
