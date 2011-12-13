% function Prob = assignmentofprodbatchestom(batches,
% capacity, costs);
%
% Creates a TOMLAB MIP problem for assignment of production batches to machines
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 18, 2005.   Last modified Oct 18, 2005.

function Prob = assignmentofprodbatchestom(batches, capacity, costs)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(batches) | isempty(capacity) | isempty(costs)
   error('One of the inputs are empty');
end

n1 = size(batches,2);  %batches.
n2 = length(capacity); %machines.
n  = n1*n2;            %machine (5) machine (5) machines (5)....until batches

% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = ones(n,1);

% Machine constr.
A1 = zeros(n1,n);
for i=1:n1
   A1(i,(i-1)*n2+1:i*n2) = ones(1,n2);
end
b_L1 = ones(n1,1);
b_U1 = b_L1;

% Batch constr.
A2 = zeros(n2,n);
for i=1:n2
   A2(i,i:n2:n-n2+i) = batches(i,:);
end
b_L2 = -inf*ones(n2,1);
b_U2 = capacity;

% Merge A, b_L, b_U
A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Objective
c   = costs(:);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production Batches to Machines', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.