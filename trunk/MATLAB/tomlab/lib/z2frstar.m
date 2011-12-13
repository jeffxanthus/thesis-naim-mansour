% z2frstar converts a table of arcs and corresponding costs in a network to
% the forward-reverse star data storage representation
%
% function [P,Z,c,T,R,u] = z2frstar(Z, C, U)
%
%	INPUT PARAMETERS
%
%   Z:	A table with arcs (i,j) node i - node j . Z is n x 2.
%       m = number of nodes = max(max(Z)). n = number of arcs.
%   C:	cost for each arc, n-vector.
%   U:	upper bounds on flow (optional)
%
%	OUTPUT PARAMETERS
%
%	P:      Pointer vector to start of each node in Z-matrix
%	Z:      Arcs outgoing from the nodes in increasing order
%	        Z(:,1) Tail. Z(:,2) Head.
%	c:      Costs related to the arcs in the Z-matrix
%	T:      Trace vector points to Z with sorting order HEAD.
%	R:      Pointer vector in T vector for each node.
%	u:	upper bounds on flow if U input or infinity

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1994-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Oct 26, 1994.  Last modified Mar 1, 2003.

function [P,Z,c,T,R,u] = z2frstar(Z, C, U)

if nargin  < 3
   U = [];
end

n = size(Z,1);
m = max(max(Z));
C = C(:); 
[y ix] = sort(Z(:,1));
Z = Z(ix,:);
c = C(ix);

if isempty(U)
   u=Inf*ones(n,1);
else
   u=U(ix);
end

% Forward star representation
P = ones(m+1,1);
k = 1;
for i = 1:n
    j = Z(i,1);
    if j > k
       while j > k
             k = k+1;
             P(k) = i;
       end
    end
end
for i = k+1:m+1
    P(i) = n+1;
end
% Reverse star representation
[y T]  = sort(Z(:,2));
R      = ones(m+1,1);
R(m+1) = n+1;
if 0
for i = 1:m
    j   = m+1-i;
    idx = find(y == j);
    if isempty(idx)
       R(j) = R(j+1);
    else
       R(j) = idx(1);
    end
end
else
k      = m;
for i=length(y):-1:1
    v = y(i);
    while k ~= v & k > 0
        R(k) = i + 1;
        k    = k - 1;
    end
end
end