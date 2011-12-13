% a2frstar converts a node-arc incidence matrix representing a network to
% the forward-reverse star data storage representation
%
% function [P,Z,c,T,R,u] = a2frstar(A, C, U)
%
%	INPUT PARAMETERS
%
%   A:	The node-arc incidence matrix. A is m x n.
%       m = number of nodes. n = number of arcs.
%   C:	cost for each arc, n-vector.
%	U:	upper bounds on flow (optional)
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
% Copyright (c) 1997-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written Nov 19, 1994.    Last modified Aug 13, 2009.

function [P,Z,c,T,R,u] = a2frstar(A, C, U)

[m n] = size(A);

if nargin < 3
   U=Inf*ones(n,1);
else
   U=U(:);
end

C = C(:); 

% Forward star representation
P = ones(m+1,1);
P(m+1) = n+1;
c = zeros(n,1);
u = zeros(n,1);
Z = zeros(n,2);
k = 1;
for i = 1:m
    idx = find(A(i,:) == 1);
    P(i) = k;
    for j = 1:length(idx)
        Z(k,:) = [i find(A(:,idx(j)) == -1)];
        c(k) = C(idx(j));
        u(k) = U(idx(j));
        k = k + 1;
    end
end
% Reverse star representation
[y T] = sort(Z(:,2));
R = ones(m+1,1);
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
k = m;
for i=length(y):-1:1
    v = y(i);
    while k ~= v & k > 0
        R(k) = i + 1;
        k    = k - 1;
    end
end
end

% MODIFICATION LOG:
%
% 001105 med 'end' added
% 090813 med  mlint check