% function Prob = routingtelephonecalls(arcs_in, arcs_out, demand_in,...
%   demand_out, demands, path_mat, path)
%
% Creates a TOMLAB MIP problem for routing telephone calls
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 22, 2005.   Last modified Nov 22, 2005.

function Prob = routingtelephonecalls(arcs_in, arcs_out, capacity, demand_in,...
   demand_out, demands, path_mat, paths)

if nargin < 7
   error('The function requires 7 inputs');
end

if isempty(arcs_in) | isempty(arcs_out) | isempty(demand_in)...
      | isempty(demand_out) | isempty(demands) | isempty(path_mat)...
      | isempty(paths) | isempty(capacity)
   error('One of the inputs are empty');
end

n1 = length(arcs_in);
n2 = length(demands);
n  = size(paths,1);   % paths

% FORMULATE PROBLEM
% No variables are binary.
IntVars   = zeros(n,1);
x_L       = zeros(n,1);
x_U       = inf*ones(n,1);

% Capacity constraint on arc
A1   = zeros(n1,n);
b_L1 = -inf*ones(n1,1);
b_U1 = capacity;
for i=1:n1
   [idx_i, idx_j] = find(path_mat(1:end,1:end-1) == arcs_in(i));
   next_nodes = path_mat(idx_i+n*idx_j);
   idx2 = find(next_nodes == arcs_out(i));
   new_idx = idx_i(idx2);   
   A1(i,new_idx) = ones(1,length(new_idx));
end

% Demand constraint
A2   = zeros(n2,n);
b_L2 = -inf*ones(n2,1);
b_U2 = demands;
for i=1:n2
   [idx_i, idx_j] = find(paths(1:end,1) == demand_in(i));
   next_nodes = paths(idx_i+n*idx_j);
   idx2 = find(next_nodes == demand_out(i));
   new_idx = idx_i(idx2);
   A2(i,new_idx) = ones(1,length(new_idx));
end

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];
c   = -ones(n,1);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Routing Telephone Calls', [], [], IntVars);

% MODIFICATION LOG
%
% 051122 med   Created.