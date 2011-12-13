% function Prob = paintproduction(prodtimes, cleantimes)
%
% Creates a TOMLAB MIP problem for paint production
%
% INPUT PARAMETERS
% prodtimes     Production times.
% cleantimes    Clean times between batches.
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 13, 2005.   Last modified Jan 2, 2006.

function Prob = paintproduction(prodtimes, cleantimes)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(cleantimes) | isempty(prodtimes)
   error('One of the inputs are empty');
end

n1 = size(cleantimes,1); % Batches, slots
n  = n1*n1 + n1;

prodtimes = prodtimes(:);

% FORMULATE PROBLEM

% All slots are integers
IntVars = [ones(n1*n1,1);zeros(n1,1)];
x_L = zeros(n,1);
x_U = [ones(n1*n1,1);inf*ones(n1,1)];

for i=1:n1+1:n1*n1
   x_U(i,1) = 0;
end

% Only one transition at a given time
A1 = zeros(n1,n);
for i=1:n1
   A1(i,(i-1)*n1+1:i*n1) = ones(1,n1);
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% Only one transition to a given batch
A2 = zeros(n1,n);
for i=1:n1
   A2(i,i:n1:n1*n1-n1+i) = ones(1,n1);
end
b_L2 = ones(n1,1);
b_U2 = ones(n1,1);

% Sub-cycle constraint
n2 = n1*(n1-1)-(n1-1);
A3 = zeros(n2,n);
counter = 1;
for i=1:n1
   for j=2:n1
      if i~=j
         A3(counter, [n1*n1+j, n1*n1+i, i+(j-1)*n1]) = [1 -1 -n1];
         counter = counter + 1;
      end
   end
end
b_L3 = (1-n1)*ones(n2,1);
b_U3 = inf*ones(n2,1);

% Add A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
c1 = repmat(prodtimes, n1, 1);
c2 = cleantimes(:);
c = [c1+c2;zeros(n1,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Paint Production', [], [], IntVars);

% MODIFICATION LOG
%
% 051010 med   Created
% 060102 med   Updated with correct model