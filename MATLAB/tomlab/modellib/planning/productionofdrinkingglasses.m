% function Prob = productionofdrinkingglasses(demand, workermax, machinemax, maxstorage...
%    , productioncost, storagecost, initialstock, finalstock, timeworker...
%    , timemachine, storagespace)
%
% Creates a TOMLAB MIP problem for production of drinking glasses
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 13, 2005.   Last modified Oct 13, 2005.

function Prob = productionofdrinkingglasses(demand, workermax, machinemax, maxstorage...
   , productioncost, storagecost, initialstock, finalstock, timeworker...
   , timemachine, storagespace)

if nargin < 11
   error('The function requires 11 inputs');
end

if isempty(demand) | isempty(workermax) | isempty(machinemax) | isempty(maxstorage)...
      | isempty(productioncost) | isempty(storagecost) | isempty(initialstock) ...
      | isempty(finalstock) | isempty(timeworker) | isempty(timemachine) ...
      | isempty(storagespace)
   error('One of the inputs are empty');
end

n1 = size(demand,1); %PRODUCTS
n2 = size(demand,2); %TIME SEGMENTS
n  = n1*n2*2; % [period-produce (6) , period-store (6) , period.... aso

% FORMULATE PROBLEM

% All slots are integers
IntVars    = zeros(n,1);
x_L        = zeros(n,1);
x_U        = inf*ones(n,1);

% Production, storage equilibrium, first period.
A1 = zeros(n1,n);

for i=1:n1
   A1(i,i:n1:2*n1-n1+i) = [1 -1];
end
b_L1 = demand(:,1)-initialstock;
b_U1 = demand(:,1)-initialstock;

% Production, storage equilibrium, all other periods.
A2 = zeros(n1*n2-n1,n);
for i=1:n1
   for j=2:n2
      A2( (i-1)*(n2-1)+(j-1) , [(j-1)*n1*2+i ; (j-1)*n1*2+n1+i; (j-2)*n1*2+n1+i] ) = [1 -1 1];
   end
end
demand(:,1) = [];
b_L2 = demand(:);
b_U2 = demand(:);

% Worker capacity constraint
A3 = zeros(n2,n);
for i=1:n2
   A3(i,(i-1)*n1*2+1:(i-1)*n1*2+n1) = timeworker';
end

b_L3 = -inf*ones(n2,1);
b_U3 = workermax*ones(n2,1);

% Machine capacity constraint
A4 = zeros(n2,n);
for i=1:n2
   A4(i,(i-1)*n1*2+1:(i-1)*n1*2+n1) = timemachine';
end

b_L4 = -inf*ones(n2,1);
b_U4 = machinemax*ones(n2,1);

% Storage constraint

A5 = zeros(n2,n);
for i=1:n2
   A5(i,(i-1)*n1*2+1+n1:(i-1)*n1*2+n1+n1) = storagespace';
end

b_L5 = -inf*ones(n2,1);
b_U5 = maxstorage*ones(n2,1);

% Final stock constraint, bounds on decision variable
x_L(n-n1+1:n,1) = finalstock;

% Add A, b_L, b_U
A = [A1;A2;A3;A4;A5];
b_L = [b_L1;b_L2;b_L3;b_L4;b_L5];
b_U = [b_U1;b_U2;b_U3;b_U4;b_U5];

% Objective

c = [repmat([productioncost;storagecost],n2,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production of Drinking Glasses', [], [], IntVars);

% MODIFICATION LOG
%
% 051017 med   Created.