% function Prob = airlinehublocation(frights, distance, hubs)
%
% Creates a TOMLAB MIP problem for airline hub location
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 7, 2005.   Last modified Dec 8, 2005.

function Prob = airlinehublocation(frights, distance, hubs)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(frights) | isempty(distance) | isempty(hubs)
   error('One of the inputs are empty');
end

n1 = size(frights,1);   % Cities
n2 = n1^4;              % i, j, k, l
n  = n2 + n1;           % i, j, k, l, then hubs
% FORMULATE PROBLEM

% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Hub constraint
A1   = [zeros(1,n2), ones(1,n1)];
b_L1 = hubs;
b_U1 = hubs;

% i and j constraint
A2   = zeros(n1^2,n);
b_L2 = ones(n1^2,1);
b_U2 = ones(n1^2,1);
for i=1:n1
   for j=1:n1
      A2((i-1)*n1+j, i+(j-1)*n1:n1^2:n2-n1^2+i+(j-1)*n1) = ones(1,n1^2);
   end
end

% Hub constraint 1
A3   = zeros(n2,n);
b_L3 = -inf*ones(n2,1);
b_U3 = zeros(n2,1);
counter = 1;
for L=1:n1
   for k=1:n1
      for j=1:n1
         for i=1:n1
            A3(counter, [counter, n2+k]) = [1 -1];
            counter = counter + 1;
         end
      end
   end
end

% Hub constraint 2
A4   = zeros(n2,n);
b_L4 = -inf*ones(n2,1);
b_U4 = zeros(n2,1);
counter = 1;
for L=1:n1
   for k=1:n1
      for j=1:n1
         for i=1:n1
            A4(counter, [counter, n2+L]) = [1 -1];
            counter = counter + 1;
         end
      end
   end
end

% Add constraints
A   = [A1;A2;A3;A4];
b_L = [b_L1;b_L2;b_L3;b_L4];
b_U = [b_U1;b_U2;b_U3;b_U4];

% Objective
c    = zeros(n,1);
counter = 1;
for L=1:n1
   for k=1:n1
      for j=1:n1
         for i=1:n1
            c(counter, 1) = frights(i,j)*(distance(i,k)+distance(L,j)+distance(k,L)*.8);
            counter = counter + 1;
         end
      end
   end
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Airline Hub Location', [], [], IntVars);

% MODIFICATION LOG
%
% 051107 med   Created.
% 051208 med   Fixed bug in hub constraint 1