% function Prob = planningtheprodofeleccomp(demand, prodcost,...
%       storagecost, initialstock, finalstock, increasecost, decreasecost);
%
% Creates a TOMLAB MIP problem for planning the production of electronic components
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 18, 2005.   Last modified Oct 18, 2005.

function Prob = planningtheprodofeleccomp(demand, prodcost,...
                storagecost, initialstock, finalstock, increasecost, decreasecost)

if nargin < 7
   error('The function requires 7 inputs');
end

if isempty(demand) | isempty(prodcost) | isempty(storagecost)...
      | isempty(initialstock) | isempty(finalstock) | isempty(increasecost)...
      | isempty(decreasecost)
   error('One of the inputs are empty');
end

n1 = size(demand,2);      % months, 6
n2 = size(demand,1);      % prods,  4 
n3 = n2*2+2;              % prods,  4, store (4), add (1), reduce (1)
n  = n1*n3;               % month (states),.. month, month, month, month, month

% FORMULATE PROBLEM

% All slots are integers
IntVars    = ones(n,1);
x_L        = [zeros(n3*(n1-1),1);zeros(n2,1);finalstock;zeros(2,1)];
% Final stock constraint
x_U        = inf*ones(n,1);

% Production equilibrium constraint at start
A1 = zeros(n2,n);
for i=1:n2
   A1(i,[i;i+n2]) = [1 -1];
end
b_L1 = demand(:,1)-initialstock;
b_U1 = b_L1;

% Production equilibrium in process
A2 = zeros(n2*(n1-1),n);
for i=1:n2
   for j=2:n1
      A2((i-1)*(n1-1)+(j-1),[n3*(j-1)+i;n3*(j-1)+n2+i;n3*(j-2)+n2+i]) = [1 -1 1];
   end
end
demand = demand(:,2:end);
demand = demand';
b_L2 = demand(:);
b_U2 = b_L2;

% Add/reduction constraint
A3 = zeros(n1-1,n);
for i=2:n1
   A3(i-1,[n3*(i-1)+1:n3*(i-1)+n2 , ...
      n3*(i-2)+1:n3*(i-2)+n2, n3*(i-1)+1+n2*2, ...
      n3*(i-1)+2+n2*2]) = [ones(1,n2), -ones(1,n2), -1, 1];
end
b_L3 = zeros(n1-1,1);
b_U3 = b_L3;

% Merge A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
mat = [prodcost;storagecost;increasecost;decreasecost];
c = repmat(mat,n1,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production of Electronic Components', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.