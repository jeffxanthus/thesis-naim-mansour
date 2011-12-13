% function Prob = depotlocation(delcosts, buildcosts, capacity, demand);
%
% Creates a TOMLAB MIP problem for depot location
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 19, 2005.   Last modified Oct 19, 2005.

function Prob = depotlocation(delcosts, buildcosts, capacity, demand, idxinf)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(delcosts) | isempty(buildcosts) | isempty(capacity)...
      | isempty(demand) | isempty(idxinf)
   error('One of the inputs are empty');
end

n1 = length(buildcosts); % DEPOTS
n2 = length(demand);     % CUSTOMERS

n  = n1+n1*n2;   % depots (exist), customer (depot1), customer (depot2)....
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = [ones(n1,1);zeros(n1*n2,1)];
x_L        = zeros(n,1);
x_U        = ones(n,1);
x_U(n1+idxinf) = 0;

% Customer constraint.
A1 = zeros(n2,n);
for i=1:n2
   A1(i,n1+i:n2:n1+n1*n2-n2+i) = ones(1,n1);
end
b_L1 = ones(n2,1);
b_U1 = ones(n2,1);

% Capacity constraint.
A2 = zeros(n1,n);
for i=1:n1
   A2(i,[i, n1+(i-1)*n2+1:n1+i*n2]) = [-capacity(i), demand'];
end
b_L2 = -inf*ones(n1,1);
b_U2 = zeros(n1,1);

% Merge A, b_L, b_U
A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Objective
delcosts = delcosts';
c   = [buildcosts;delcosts(:)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Depot Location', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.