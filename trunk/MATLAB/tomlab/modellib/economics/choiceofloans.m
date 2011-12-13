% function Prob = choiceofloans(rates, costs, maxloan, loanlength)
%
% Creates a TOMLAB MIP problem for location of gsm transmitters
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = choiceofloans(rates, costs, maxloan, loanlength)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(rates) | isempty(costs) | isempty(maxloan) | isempty(loanlength)
   error('One of the inputs are empty');
end

n1 = size(rates,1);       %Banks
n2 = size(rates,2);       %Shops
n  = n1 * n2; % Banks (shop 1), Banks (shop 2), Banks (shop 3)

% FORMULATE PROBLEM
% No variables are binary.
IntVars   = zeros(n,1);
x_L       = zeros(n,1);
x_U       = inf*ones(n,1);

% Shop constraints
A1   = zeros(n2,1);
b_L1 = costs(:);
b_U1 = costs(:);
for i=1:n2
   A1(i,i:n1:n-n1+i) = ones(1,n1);
end

% Bank constraints
A2   = zeros(n1,1);
b_L2 = -inf*ones(n1,1);
b_U2 = maxloan*ones(n1,1);
for i=1:n1
   A2(i,(i-1)*n2+1:i*n2) = ones(1,n2);
end

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];
rates = rates';
rates = rates(:);

c   = (rates/100)./(1-(1+rates/100).^(-loanlength));
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Choice of Loans', [], [], IntVars);

% MODIFICATION LOG
%
% 051201 med   Created.