% function Prob = carrental(cost, distance, stock, demand);
%
% Creates a TOMLAB MIP problem for car rental
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 19, 2005.   Last modified Oct 19, 2005.

function Prob = carrental(cost, distance, stock, demand)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(cost) | isempty(distance) | isempty(stock) | isempty(demand)
   error('One of the inputs are empty');
end

idx_excess = find(stock-demand > 0);
n_excess   = length(idx_excess);

idx_need   = find(stock-demand < 0);
n_need     = length(idx_need);

n  = n_excess*n_need;  % from 1 to others are first column
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = inf*ones(n,1);

% Excess constraint
A1 = zeros(n_excess,n);
for i=1:n_excess
   A1(i,i:n_excess:n-n_excess+i) = ones(1,n_need);
end
b_L1 = stock(idx_excess) - demand(idx_excess);
b_U1 = b_L1;

% Need constraint
A2 = zeros(n_need,n);
for i=1:n_need
   A2(i,(i-1)*n_excess+1:i*n_excess ) = ones(1,n_excess);
end
b_L2 = demand(idx_need) - stock(idx_need);
b_U2 = b_L2;

% Merge A, b_L, b_U
A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Objective
c   = cost*distance(idx_excess,idx_need);
c = c(:);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Car Rental', [], [], IntVars);
Prob.user.idx_excess = idx_excess;
Prob.user.idx_need = idx_need;

% MODIFICATION LOG
%
% 051018 med   Created.