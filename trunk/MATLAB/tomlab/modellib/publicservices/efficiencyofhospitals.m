% function Prob = efficiencyofhospitals(resources, services, hospidx)
%
% Creates a TOMLAB MILP problem for efficiency of hospitals
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 6, 2005.   Last modified Dec 6, 2005.

function Prob = efficiencyofhospitals(resources, services, hospidx)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(resources) | isempty(services) | isempty(hospidx)
   error('One of the inputs are empty');
end

n1   = size(resources,2);    %Coef
n2   = size(services,1);     %Services
n3   = size(resources,1);    %Resources
n4   = 1;                    %Eff

n    = n1+n2+n3+n4; %Coef, Servs, Res, Eff

% FORMULATE PROBLEM
% No variables are binary
IntVars   = zeros(n,1);
x_L       = zeros(n,1);
x_U       = inf*ones(n,1);

% Coef constraints
b_L1  = 1;
b_U1  = 1;
A1    = [ones(1,n1), zeros(1,n-n1)];

% Service constraint
b_L2  = zeros(n2,1);
b_U2  = zeros(n2,1);
A2    = zeros(n2,n);
for i=1:n2
   A2(i,[1:n1, n1+i]) = [services(i,:), -1];
end

% Resource constraint
b_L3  = zeros(n3,1);
b_U3  = zeros(n3,1);
A3    = zeros(n3,n);
for i=1:n3
   A3(i,[1:n1,n1+n2+i]) = [resources(i,:), -1];
end

% Service indicators, greater than ficticious ones
x_L(n1+1:n1+n2,1) = services(:,hospidx);

% Efficiency relationship
b_L4  = -inf*ones(n3,1);
b_U4  = zeros(n3,1);
A4    = zeros(n3,n);
for i=1:n3
   A4(i,[n1+n2+i,n]) = [1,-sum(resources(i,hospidx))];
end

% Merge constraints
A   = [A1;A2;A3;A4];
b_L = [b_L1;b_L2;b_L3;b_L4];
b_U = [b_U1;b_U2;b_U3;b_U4];

c = zeros(n,1);
c(end,1) = 1;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Efficiency of Hospitals', [], [], IntVars);

% MODIFICATION LOG
%
% 051206 med   Created.