% function Prob = choosingthemodeoftransport(arcs_out, arcs_in, arcs_min,
% arcs_max, arcs_cost);
%
% Creates a TOMLAB MIP problem for choosing the mode of transport
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 19, 2005.   Last modified Oct 19, 2005.

function Prob = choosingthemodeoftransport(arcs_out, arcs_in, arcs_min, arcs_max, arcs_cost, minflow)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(arcs_out) | isempty(arcs_in) | isempty(arcs_min)...
      | isempty(arcs_max) | isempty(arcs_cost) | isempty(minflow)
   error('One of the inputs are empty');
end

n  = length(arcs_out); % ARCS
n1 = length(unique([arcs_in ; arcs_out]));
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = arcs_min;
x_U        = arcs_max;

% Node constraint, except for source and sink

A1 = zeros(n1-2,n);
for i=2:n1-1
   idx_in  = find(arcs_in  == i);
   idx_out = find(arcs_out == i);
   A1(i-1,[idx_in;idx_out]) = [ones(1,length(idx_in)), -ones(1,length(idx_out))];
end
b_L1 = zeros(n1-2,1);
b_U1 = b_L1;

% Source constraint
A2 = zeros(1,n);
idx_in  = find(arcs_out  == 1); %1 is always source
A2(1,idx_in) = ones(1,length(idx_in));

b_L2 = minflow;
b_U2 = b_L2;

% Merge A, b_L, b_U
A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Objective
c   = arcs_cost;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Choosing the Mode of Transport', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.