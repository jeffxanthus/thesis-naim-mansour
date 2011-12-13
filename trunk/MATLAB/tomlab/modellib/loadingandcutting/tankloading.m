% function Prob = tankloading(capacity, products)
%
% Creates a TOMLAB MIP problem for tank loading.
%
% INPUT PARAMETERS
% capacity     Capacity of each tank
% products     Amount to fill of each product
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = tankloading(capacity, products)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(capacity) | isempty(products)
   error('One of the inputs are empty');
end

capacity = capacity(:);
products = products(:);

% FORMULATE PROBLEM
n1 = length(capacity);
n2 = length(products);
n = n1*n2;

c = repmat(capacity,n2,1); % - since minimizing
x_L = zeros(n,1);
x_U = ones(n,1);
IntVars = ones(n,1);

% Load constraint, the liquids have to go
A1 = zeros(n2,n);
for i=1:n2
   A1(i,(i-1)*n1+1:i*n1) = capacity';
end
b_L1 = products;
b_U1 = inf*ones(n2,1);

% Only one product per tank
A2 = zeros(n1,n);
for i=1:n1
   A2(i,i:n1:n-n1) = 1;
end
b_L2 = -inf*ones(n1,1);
b_U2 = ones(n1,1);

% Add constraints
A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Tank Loading',...
   [],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.