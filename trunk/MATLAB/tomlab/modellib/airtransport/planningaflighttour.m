% function Prob = planningaflighttour(distance)
%
% Creates a TOMLAB MIP problem for planning a flight tour
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 7, 2005.   Last modified Nov 7, 2005.

function Prob = planningaflighttour(distance)

if nargin < 1
   error('The function requires 1 input');
end

if isempty(distance)
   error('One of the inputs are empty');
end

n1 = size(distance,1);  % Cities
n  = n1 * n1;           % i and j

% FORMULATE PROBLEM
% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Cycle constraints
A1   = zeros(n1,n);
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);
for i=1:n1
   for j=1:n1
      if i~=j
         A1(i,(j-1)*n1+i) = 1;
      end
   end
end
A2   = zeros(n1,n);
b_L2 = ones(n1,1);
b_U2 = ones(n1,1);
for i=1:n1
   for j=1:n1
      if i~=j
         A2(i,(i-1)*n1+j) = 1;
      end
   end
end

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

c    = distance(:);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Planning a Flight Tour', [], [], IntVars);

% MODIFICATION LOG
%
% 051107 med   Created.