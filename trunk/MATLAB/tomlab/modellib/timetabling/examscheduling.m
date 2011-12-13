% function Prob = examscheduling(incompatmat, slots)
%
% Creates a TOMLAB MILP problem for exam scheduling
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = examscheduling(incompatmat, slots)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(incompatmat) | isempty(slots)
   error('One of the inputs are empty');
end

n1    = size(incompatmat, 1);    %exams
n2    = slots;                   %8 time slots

n     = n1*n2;        %exams, exams, exams...

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Exam constraint
b_L1 = ones(n1,1);
b_U1 = b_L1;
A1   = zeros(n1,n);
for i=1:n1
   A1(i,i:n1:n-n1+i) = ones(1,n2);
end

% Incompatibility constr.
ncon = nnz(triu(incompatmat));
b_L2 = -inf*ones(ncon*n2,1);
b_U2 = ones(ncon*n2,1);
A2   = zeros(ncon*n2,n);
count = 1;
for i=1:n1-1
   for j=i+1:n1
      if incompatmat(i,j) == 1
         for k=1:n2
            A2(count,[(k-1)*n1+i,(k-1)*n1+j]) = ones(1,2);
            count = count + 1;
         end
      end
   end
end

% Merge of constraints
A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

c    = zeros(n,1);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Exam Scheduling', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.