% function Prob = jobshopscheduling(proctimes, flow, final, bigM)
%
% Creates a TOMLAB MIP problem for job shop scheduling
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Jan 2, 2006.   Last modified Jan 2, 2006.

function Prob = jobshopscheduling(proctimes, flow, final, bigM)

if nargin < 4
   error('The function requires 2 inputs');
end

if isempty(proctimes) | isempty(flow) | isempty(final) | isempty(bigM)
   error('One of the inputs are empty');
end

n1 = size(proctimes,2); % Number of jobs     (j in JOBS)
n2 = size(proctimes,1); % Number of machines (m in MACH)
n3 = 1;                 % Finish time        (scalar)

n = n1*n2+n3;

% machine 1(jobs), machine 2(jobs).... finish

% FORMULATE PROBLEM

% Disjunctive constraint
numconstr = zeros(n2,1);
numvbls   = zeros(n2,1);
for i=1:n2
   [idx1, idx2] = find(flow == i);
   numvbls(i,1) = length(idx1);
   numconstr(i,1) = prod(1:length(idx1));
end

% Expand problem
nextra = sum(numconstr);
n = n + nextra;

% All slots are integers
IntVars = ones(n,1);
x_L = zeros(n,1);
x_U = [inf*ones(n-nextra,1);ones(nextra,1)];

% Finish constraint
n4 = length(final);
A1 = zeros(n4,n);
b_L1 = zeros(n4,1);
b_U1 = inf*ones(n4,1);

for i=1:n4
   A1(i,[n-nextra, (final(i)-1)*n1+i]) = [1 -1];
   b_L1(i,1) = proctimes(final(i), i);
end

% Flow constraints.
n5 = nnz(proctimes) - length(final);
A2 = zeros(n5,n);
b_L2 = -inf*ones(n5,1);
b_U2 = zeros(n5,1);
counter = 1;
for j=2:n1
   for m=1:n2
      if ~(flow(m,j) == 0)
         A2(counter,[(flow(m,j-1)-1)*n1+m, (flow(m,j)-1)*n1+m]) = [1,-1];
         b_U2(counter, 1) = -proctimes(flow(m,j-1),m);
         counter = counter + 1;
      end
   end
end

% Disjunctive constraints, need good way to add.
vec11 = [1 2 1 3 2 3]';
vec12 = [2 1 3 1 3 2]';
vec13 = [1 1 2 2 3 3]';
vec21 = [2 3]';
vec22 = [3 2]';
vec23 = [4 4]';
vec31 = [1 2 1 3 2 3]';
vec32 = [2 1 3 1 3 2]';
vec33 = [5 5 6 6 7 7]';

n6 = sum(numconstr);
b_L3 = -inf*ones(n6, 1);
A3 = zeros(1,n);
for i=1:numconstr(1)
   if mod(i,2)
      A3(i, [vec11(i), vec12(i), n-nextra+vec13(i)]) = [1 -1 bigM];
      b_U3(i,1) = bigM - proctimes(1,vec11(i));
   else
      A3(i, [vec11(i), vec12(i), n-nextra+vec13(i)]) = [1 -1 -bigM];
      b_U3(i,1) = -proctimes(1,vec11(i));
   end
end

for i=1:numconstr(2)
   if mod(i,2)
      A3(end+1, 3+[vec21(i), vec22(i), n-nextra+vec23(i)]) = [1 -1 bigM];
      b_U3(end+1,1) = bigM - proctimes(1,vec21(i));
   else
      A3(end+1, 3+[vec21(i), vec22(i), n-nextra+vec23(i)]) = [1 -1 -bigM];
      b_U3(end+1,1) = -proctimes(1,vec21(i));
   end
end

for i=1:numconstr(3)
   if mod(i,2)
      A3(end+1, 6+[vec31(i), vec32(i), n-nextra+vec33(i)]) = [1 -1 bigM];
      b_U3(end+1,1) = bigM - proctimes(1,vec31(i));
   else
      A3(end+1, 6+[vec31(i), vec32(i), n-nextra+vec33(i)]) = [1 -1 -bigM];
      b_U3(end+1,1) = -proctimes(1,vec31(i));
   end
end

% Add A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
c = zeros(n,1);
c(n1*n2+n3,1) = 1;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Job Shop Scheduling', [], [], IntVars);

% MODIFICATION LOG
%
% 060102 med   Created.