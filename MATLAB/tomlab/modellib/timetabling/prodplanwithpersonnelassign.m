% function Prob = prodplanwithpersonnelassig(profit, capacity, timemat,...
%   transfermat, maxtransfer)
%
% Creates a TOMLAB MILP problem for production planning with personnel assignment
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = prodplanwithpersonnelassign(profit, capacity, timemat,...
   transfermat, maxtransfer, flag)

if nargin < 6
   error('The function requires 6 inputs');
end

if isempty(profit) | isempty(capacity) | isempty(timemat) | isempty(transfermat)...
      | isempty(maxtransfer) | isempty(flag)
   error('One of the inputs are empty');
end

if flag == 0
   n1    = length(profit);    %products
   n2    = length(capacity);  %production lines
   n     = n1;        %amount

   % FORMULATE PROBLEM
   % No variables are binary
   IntVars   = zeros(n,1);
   x_L       = zeros(n,1);
   x_U       = inf*ones(n,1);

   % Capacity constraint
   b_L = -inf*ones(n2,1);
   b_U = capacity;
   A   = zeros(n2,n);
   for i=1:n2
      A(i,1:n1) = timemat(:,i)';
   end
   c    = -profit;
end

if flag == 1
   n1    = length(profit);    %products
   n2    = length(capacity);  %production lines
   n     = n1 + n2*n2 + n2;   %amount, transfers, hours

   % FORMULATE PROBLEM
   % No variables are binary
   IntVars   = zeros(n,1);
   x_L       = zeros(n,1);
   x_U       = inf*ones(n,1);

   % Line constraint
   b_L1 = -inf*ones(n2,1);
   b_U1 = zeros(n2,1);
   A1  = zeros(n2,n);
   for i=1:n2
      A1(i,[1:n1, n1+n2*n2+i]) = [timemat(:,i)',-1];
   end
   
   % Line transfer constraint
   b_L2 = -capacity;
   b_U2 = -capacity;
   A2  = zeros(n2,n);
   for i=1:n2
      val1 = transfermat(:,i);
      idx1 = find(val1 == 1);
      idx1 = idx1 + (i-1)*n2+n1;
      
      val2 = transfermat(i,:)';
      idx2 = find(val2 == 1);
      idx2 = (idx2-1)*n2+i+n1;
      A2(i,[idx1', n1+n2*n2+i, idx2']) = [ones(1,length(idx1)),-1,-ones(1,length(idx2))];
   end
   
   % Line max transfer constraint
   b_L3 = -inf*ones(n2,1);
   b_U3 = maxtransfer;
   A3  = zeros(n2,n);
   for i=1:n2
      val2 = transfermat(i,:)';
      idx2 = find(val2 == 1);
      idx2 = (idx2-1)*n2+i+n1;
      A3(i,idx2) = ones(1,length(idx2));
   end
   
   % Merge constraints
   A   = [A1;A2;A3];
   b_L = [b_L1;b_L2;b_L3];
   b_U = [b_U1;b_U2;b_U3];
   c    = [-profit;zeros(n2*n2+n2,1)];
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production Planning with Personnel Assignment', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.