% function Prob = planningtheprodoffiberglass(capacity, demand,...
%       prodcost, storcost);
%
% Creates a TOMLAB MIP problem for planning the production of fiber glass
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 18, 2005.   Last modified Oct 18, 2005.

function Prob = planningtheprodoffiberglass(capacity, demand,...
       prodcost, storcost)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(capacity) | isempty(demand) | isempty(prodcost)...
      | isempty(storcost)
   error('One of the inputs are empty');
end

n1 = size(demand,1);
n  = n1*2-1;

% FORMULATE PROBLEM

% No variables are integers (cost flow system)
IntVars    = zeros(n,1);
x_L        = zeros(n,1);

% Flow constraints
vec        = [capacity, capacity];
vec        = vec';
vec        = vec(:);
x_U        = vec(1:end-1,1);

% First node constraint
A1 = zeros(1,n);
A1(1,[1;2]) = [1 -1];
b_L1 = demand(1);
b_U1 = b_L1;

% Constraints for all other nodes, except final
A2 = zeros(n1-2,n);
for i=1:n1-2
   A2(i,[i*2, i*2+1, i*2+2]) = [1 1 -1];
end
b_L2 = demand(2:end-1,1);
b_U2 = b_L2;

% First node constraint
A3 = zeros(1,n);
A3(1,[n-1;n]) = [1 1];
b_L3 = demand(end,1);
b_U3 = b_L3;

% Merge A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
vec        = [prodcost, storcost];
vec        = vec';
vec        = vec(:);
c          = vec(1:end-1,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production of Fiber Glass', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.