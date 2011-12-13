% function Prob = schedulingflightlandings(costs, mintimes)
%
% Creates a TOMLAB MIP problem for scheduling flight landings
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 6, 2005.   Last modified Nov 6, 2005.

function Prob = schedulingflightlandings(costs, mintimes)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(costs) | isempty(mintimes)
   error('One of the inputs are empty');
end

n1 = size(costs,2);     % Planes, land(p) time when landing
n2 = n1*n1;             % Binary vector, landing p precedes q
n  = n1 + n2 + n1 + n1; % land(p), bin vec, early(p), late(p)
ncon = sum(1:n1-1);
BIG = 1e4;
% FORMULATE PROBLEM

% All variables are integers.
IntVars   = ones(n,1);
x_L       = [costs(1,:)';zeros(n2+2*n1,1)];
x_U       = [costs(3,:)';ones(n2,1);costs(2,:)'-costs(1,:)';costs(3,:)'-costs(2,:)'];

for q=1:n1-1            % q is less than p
   for p=q+1:n1
      x_U(n1+(p-1)*n1+q,1) = 0;
   end
end

for i=1:n1+1:n2         % Variables in the diagonal
   x_U(n1+i,1) = 0;
end

% Disjunctive constr 1
A1   = zeros(ncon,n);
b_L1 = -inf*ones(ncon,1);
b_U1 = zeros(ncon,1);
counter = 1;
for q=1:n1-1            % q is less than p
   for p=q+1:n1
      A1(counter,[p, q, n1+(q-1)*n1+p]) = [1 -1 -BIG];
      b_U1(counter, 1) = -mintimes(p,q);
      counter = counter + 1;
   end
end

% Disjunctive constr 2
A2   = zeros(ncon,n);
b_L2 = -inf*ones(ncon,1);
b_U2 = zeros(ncon,1);
counter = 1;
for p=1:n1-1            % p is less than q
   for q=p+1:n1 
      A2(counter,[p, q, n1+(p-1)*n1+q]) = [1 -1 BIG];
      b_U2(counter, 1) = BIG-mintimes(p,q);
      counter = counter + 1;
   end
end

% Landing constraint
A3   = zeros(n1,n);
b_L3 = costs(2,:)';
b_U3 = costs(2,:)';
for i=1:n1
   A3(i,[i, n1+n2+i, n1+n2+n1+i]) = [1 1 -1];
end

% Add constraints
A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
c                     = zeros(n,1);
c(n1+n2+1:n1+n2+n1,1) = costs(4,:)';
c(n1+n2+n1+1:n,1)     = costs(5,:)';

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Scheduling Flight Landings', [], [], IntVars);

% MODIFICATION LOG
%
% 051106 med   Created.