% function Prob = constructionofacablednetwork(distances)
%
% Creates a TOMLAB MIP problem for construction of a cabled network
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 22, 2005.   Last modified Nov 22, 2005.

function Prob = constructionofacablednetwork(distances)

if nargin < 1
   error('The function requires 1 inputs');
end

if isempty(distances)
   error('One of the inputs are empty');
end

n1 = size(distances, 1);
n2 = n1-1;
n3 = sum(1:n2);  % binary vbls
n  = n3 + n1;

% FORMULATE PROBLEM
% Some variables are binary.
IntVars   = [ones(n3,1);zeros(n1,1)];
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Terminal connections constraint
A1   = ones(1,n);
b_L1 = n2;
b_U1 = n2;

% At least one terminal connected to other
A2   = zeros(n2,n);
b_L2 = ones(n2,1);
b_U2 = inf*ones(n2,1);
counter = 0;
idx_anti_in = [];
idx_anti_out = [];
for i=n2:-1:1
   idx1 = sum(i:n2)-n2+1+counter;
   idx2 = idx1 + i-1;
   idx  = idx1:idx2;
   A2(i,idx) = ones(1,length(idx));
   counter = counter + 1;
   idx_anti_in = [idx_anti_in;counter*ones(length(idx),1)];
   idx_anti_out = [idx_anti_out;[counter+1:n1]'];
end

% Anti-cycling constraint
A3   = zeros(n3,n);
b_L3 = -inf*ones(n3,1);
b_U3 = n2*ones(n3,1);
for i=1:n3
   idx = [i; n3+idx_anti_in(i); n3+idx_anti_out(i)];
   A3(i,idx) = [n2, 1, -1];
end

A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

distances = distances(tril(distances) ~= 0);

c   = [distances(:);zeros(n1,1)];
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Construction of a Cabled Network', [], [], IntVars);

% MODIFICATION LOG
%
% 051122 med   Created.