% tp2np Convert Transportation LP to Network Programming formulation
%       using supply, demand and distance matrix
%
% function Net = tp2np(s, d, C)
%
%                    m = M + N; n = M * N
%
%   Input:
%           s        Supply, N by 1
%           d        Demand, M by 1
%           dist     Distance matrix, M by N
%   Output:
%      Net  Structure storing network in Forward-Reverse Star Representation
%
%           Fields used:
%
%		P:      Pointer vector to start of each node in Z-matrix
%		Z:      Arcs outgoing from the nodes in increasing order
%		        Z(:,1) Tail. Z(:,2) Head.
%		c:      Costs related to the arcs in the Z-matrix
%		T:      Trace vector points to Z with sorting order HEAD.
%		R:      Pointer vector in T vector for each node.
%		u:	upper bounds on flow if U input or infinity

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Mar 1, 2003.   Last modified Mar 1, 2003.

function Net = tp2np(s, d, C)

m  = length(s);
n  = length(d);
mn = m*n;
C = C';
c = C(:);

Z = zeros(mn,2);
k = 0;
for i = 1:m
    for j=1:n
        k      = k + 1;
        Z(k,1) = i;
        Z(k,2) = m+j;
    end
end

U = [];

[Net.P,Net.Z,Net.c,Net.T,Net.R,Net.u] = z2frstar(Z, c, U);

Net.b = [s(:);-d(:)];

% MODIFICATION LOG
%
% 030301 hkh Written
