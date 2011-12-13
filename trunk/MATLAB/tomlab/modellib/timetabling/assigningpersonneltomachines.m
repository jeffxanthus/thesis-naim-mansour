% function Prob = assigningpersonneltomachines(prodmat, flag)
%
% Creates a TOMLAB MILP problem for assigning personnel to machines
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 2, 2005.   Last modified Dec 2, 2005.

function Prob = assigningpersonneltomachines(prodmat, flag)

if nargin < 2
   error('The function requires 2 input');
end

if isempty(prodmat)
   error('One of the inputs are empty');
end

if flag == 0

   n1   = size(prodmat,1);  %workers
   n2   = size(prodmat,2);  %machines
   n    = n1*n2;            %workers(mach 1)

   % FORMULATE PROBLEM
   % All variables are binary
   IntVars   = ones(n,1);
   x_L       = zeros(n,1);
   x_U       = ones(n,1);

   % Worker constraints
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);
   A1   = zeros(n1,n);
   for i=1:n1
      A1(i,i:n2:n-n1+i) = ones(1,n2);
   end

   % Machine constraints
   b_L2 = ones(n2,1);
   b_U2 = ones(n2,1);
   A2   = zeros(n2,n);
   for i=1:n2
      A2(i,(i-1)*n1+1:i*n1) = ones(1,n1);
   end

   A   = [A1;A2];
   b_L = [b_L1;b_L2];
   b_U = [b_U1;b_U2];

   c   = -prodmat(:);

end

if flag == 1

   n1   = size(prodmat,1);  %workers
   n2   = size(prodmat,2);  %machines
   n3   = 1; %pmin
   n    = n1*n2+n3;            %workers(mach 1)

   % FORMULATE PROBLEM
   % All variables are binary
   IntVars   = [ones(n-n3,1);0];
   x_L       = zeros(n,1);
   x_U       = ones(n,1); x_U(end,1) = inf;

   % Worker constraints
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);
   A1   = zeros(n1,n);
   for i=1:n1
      A1(i,i:n2:n-n3-n1+i) = ones(1,n2);
   end

   % Machine constraints
   b_L2 = ones(n2,1);
   b_U2 = ones(n2,1);
   A2   = zeros(n2,n);
   for i=1:n2
      A2(i,(i-1)*n1+1:i*n1) = ones(1,n1);
   end
   
   %Productivity bounds 1
   b_L3 = zeros(n1,1);
   b_U3 = inf*ones(n1,1);
   A3   = zeros(n1,n);
   for i=1:n1
      A3(i, [(i-1)*n2+1:i*n2, n]) = [prodmat(:,i)', -1];
   end

   A   = [A1;A2;A3];
   b_L = [b_L1;b_L2;b_L3];
   b_U = [b_U1;b_U2;b_U3];

   c   = [zeros(n-n3,1);-1];

end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Assigning Personnel to Machines', [], [], IntVars);

% MODIFICATION LOG
%
% 051202 med   Created.