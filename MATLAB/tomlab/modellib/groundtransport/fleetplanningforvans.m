% function Prob = fleetplanningforvans(demand, initialsupply,
% contractlength, contractcost);
%
% Creates a TOMLAB MIP problem for fleet planning for vans
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 21, 2005.   Last modified Oct 21, 2005.

function Prob = fleetplanningforvans(demand, initialsupply, contractlength, contractcost)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(demand) | isempty(initialsupply) | isempty(contractlength)...
      | isempty(contractcost)
   error('One of the inputs are empty');
end

n1 = length(contractlength);   % contracts
n2 = length(demand);           % months

n  = n1*n2;              % month (contracts 1), month (mode 2)...

% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = inf*ones(n,1);

% Demand constraint
A = zeros(n2,n);
for i=1:n2
   for j=1:n1
      low = max(1,i-contractlength(j)+1);
      high = min(i,n2-contractlength(j)+1);
      for k=low:high
         A(i,(k-1)*n1+j) = 1;
      end
   end
end
b_L = demand-initialsupply;
b_U = inf*ones(n2,1);

% Objective
c   = repmat(contractcost,n2,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Fleet Planning for Vans', [], [], IntVars);

% MODIFICATION LOG
%
% 051021 med   Created.