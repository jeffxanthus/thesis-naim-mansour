% function Prob = riggingelections(districts, votes, population,
% districtsnum)
%
% Creates a TOMLAB MILP problem for rigging elections
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = riggingelections(districts, votes, population, districtsnum)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(districts) | isempty(votes) | isempty(population) | isempty(districtsnum)
   error('One of the inputs are empty');
end

n    = size(districts, 1);  %possible districts
n1   = length(votes);
% Calculate majority
d = zeros(n,1);
for i=1:n
   idx = find(districts(i,:) == 1);
   if sum(votes(idx))/sum(population(idx)) > 0.5
      d(i,1) = 1;
   end
end

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Districts chosen
b_L1  = districtsnum;
b_U1  = districtsnum;
A1    = ones(1,n);

% Quarter only once
b_L2  = ones(n1,1);
b_U2  = b_L2;
A2    = ones(n1,n);
for i=1:n1
   A2(i,:) = districts(:,i)';
end

% Merge constraints
A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];
c = -d; 

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Rigging Elections', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.