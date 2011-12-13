% function Prob = publicitycampaign(budget, people, costs, maxuse, quality,
%     minpeople)
%
% Creates a TOMLAB MIP problem for publicity campaign
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = publicitycampaign(budget, people, costs, maxuse, quality, minpeople)

if nargin < 6
   error('The function requires 6 inputs');
end

if isempty(budget) | isempty(people) | isempty(costs) | isempty(maxuse) ...
      | isempty(quality) | isempty(minpeople)
   error('One of the inputs are empty');
end

n  = length(quality); % use m

% FORMULATE PROBLEM
% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = maxuse(:);

% Budget constraints
A1   = costs';
b_L1 = -inf;
b_U1 = budget;

% Reach constraint
A2   = people';
b_L2 = minpeople;
b_U2 = inf;

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

c   = -quality;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Publicity Campaign', [], [], IntVars);

% MODIFICATION LOG
%
% 051201 med   Created.