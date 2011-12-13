% function Prob = backingupfiles(maxuse, capacity, sizes)
%
% Creates a TOMLAB MIP problem for backing up files
%
% INPUT PARAMETERS
% capacity     Capacity of each object (one size)
% sizes        Sizes of objects to put.
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = backingupfiles(maxuse, capacity, sizes)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(maxuse) | isempty(capacity) | isempty(sizes)
   error('One of the inputs are empty');
end

maxuse   = maxuse(:);
capacity = capacity(:);
sizes    = sizes(:);

% FORMULATE PROBLEM
n1 = maxuse;
n2 = length(sizes);
n3 = n1*n2;
n = n3 + n1; % which unit to save on, plus units used

c = [zeros(n3,1);ones(n1,1)]; %
x_L = zeros(n,1);
x_U = ones(n,1);
IntVars = ones(n,1);

% An item must be stored on one unit
A1 = zeros(n2,n);
for i=1:n2
   A1(i,i:n2:n3-n2+i) = 1;
end
b_L1 = ones(n2,1);
b_U1 = ones(n2,1);

% Stay below capacity of storage unit
A2 = zeros(n1,n);
for i=1:n1
   A2(i,(i-1)*n2+1:i*n2) = sizes';
   A2(i,n3+i) = -capacity;
end
b_L2 = -inf*ones(n1,1);
b_U2 = zeros(n1,1);

% Add constraints
A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Backing up files',...
   [],[],IntVars);

% MODIFICATION LOG
% 051007 med   Created.