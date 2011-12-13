% function Prob = networkreliability(distance)
%
% Creates a TOMLAB MIP problem for network reliability
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 7, 2005.   Last modified Nov 7, 2005.

function Prob = networkreliability(arcs_in, arcs_out)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(arcs_in) | isempty(arcs_out)
   error('One of the inputs are empty');
end

n1 = length(arcs_in);   % arcs_in
n  = 2*n1;              % bi-directional
nodes = length(unique([arcs_in; arcs_out]));
arcs_in_new = [arcs_in;arcs_out];
arcs_out = [arcs_out;arcs_in];
arcs_in  = arcs_in_new;

% FORMULATE PROBLEM
% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Kirchhoff's law, 1 source and 1 sink
A1   = zeros(nodes-2,n);
b_L1 = zeros(nodes-2,1);
b_U1 = zeros(nodes-2,1);
for i=1:nodes-2 %Assuming 10 and 11 are source and sink
   idx1 = find(arcs_in == i);
   idx3 = find(arcs_out == i);
   A1(i,[idx1;idx3]) = [ones(1, length(idx1)), -ones(1,length(idx3))];
end

% Limit flow in intermediate node
A2   = zeros(nodes-2,n);
b_L2 = -inf*ones(nodes-2,1);
b_U2 = ones(nodes-2,1);
for i=1:nodes-2 %Assuming 10 and 11 are source and sink
   idx1 = find(arcs_in == i);
   A2(i,idx1) = ones(1, length(idx1));
end

% Source has no incoming flows
A3   = zeros(1,n);
b_L3 = 0;
b_U3 = 0;
idx1 = find(arcs_out == nodes-1);
A3(1,idx1) = 1;

A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

idx1 = find(arcs_in == (nodes-1));

c    = zeros(n,1);
c(idx1,1) = -1;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Network Reliability', [], [], IntVars);

% MODIFICATION LOG
%
% 051107 med   Created.