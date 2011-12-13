% function Prob = planningtheprodofbicycles(normcapacity, extracapacity, normcost,...
%    extracost, demand, startstock, storagecost)
%
% Creates a TOMLAB MIP problem for planning the production of bicycles
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 13, 2005.   Last modified Oct 13, 2005.

function Prob = planningtheprodofbicycles(normcapacity, extracapacity, normcost,...
   extracost, demand, startstock, storagecost)

if nargin < 7
   error('The function requires 7 inputs');
end

if isempty(normcapacity) | isempty(extracapacity) | isempty(normcost) | isempty(extracost)...
      | isempty(demand) | isempty(startstock) | isempty(storagecost)
   error('One of the inputs are empty');
end

n1 = length(demand);
n  = n1*3; % pnorm, pover, store

% FORMULATE PROBLEM

% All slots are integers
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = [normcapacity*ones(n1,1) ; extracapacity*ones(n1,1) ; inf*ones(n1,1)];

% Every process has to go to one machine
A = zeros(n1,n);
A(1,1:n1:n-n1+1) = [1 1 -1];
for i=2:n1
   A(i,[[i:n1:n-n1+i]';n-n1-1+i]) = [1 1 -1 1];
end
b_L = demand;
b_U = demand;
b_L(1,1) = b_L(1,1) - startstock;
b_U(1,1) = b_U(1,1) - startstock;

% Objective

c = [normcost*ones(n1,1);extracost*ones(n1,1);storagecost*ones(n1,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Planning the Production of Bicycles', [], [], IntVars);

% MODIFICATION LOG
%
% 051017 med   Created.