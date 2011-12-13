% function Prob = locationofgsmtransmitters(budget, cost, connections, ...
%   population);
%
% Creates a TOMLAB MIP problem for location of gsm transmitters
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 30, 2005.   Last modified Nov 30, 2005.

function Prob = locationofgsmtransmitters(budget, cost, connections, ...
   population)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(budget) | isempty(cost) | isempty(connections) | isempty(population)
   error('One of the inputs are empty');
end

n1 = length(cost);        %number of locations to build
n2 = length(population);  %areas to cover
n  = n1 + n2; % locations, areas

% FORMULATE PROBLEM
% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Budget constraint
A1   = [cost', zeros(1,n2)];
b_L1 = -inf;
b_U1 = budget;

% Coverage constraint
A2 = [connections', -eye(n2)];
b_L2 = zeros(n2,1);
b_U2 = inf*ones(n2,1);

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

c   = [zeros(n1,1);-population(:)];
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Location of GSM Transmitters', [], [], IntVars);

% MODIFICATION LOG
%
% 051130 med   Created.