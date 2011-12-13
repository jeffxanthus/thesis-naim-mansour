% function Prob = flowshopscheduling(nummachines, proctimes);
%
% Creates a TOMLAB MIP problem for flow shop scheduling
%
% INPUT PARAMETERS
% nummachines   Number of machines for each task.
% proctimes     Processing time for each task (column wise for each unit).
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 13, 2005.   Last modified Jan 2, 2006.

function Prob = flowshopscheduling(nummachines, proctimes)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(nummachines) | isempty(proctimes)
   error('One of the inputs are empty');
end

nummachines  = nummachines(:);

n1 = size(proctimes,2); % Number of jobs     (j in JOBS)
n2 = size(proctimes,2); % Number of ranks    (k in RANKS)
n3 = nummachines;       % Number of machines (m in MACH)

n = n1*n2 + n3*(n1-1) + (n3-1)*n2;

% Variables, rank(j,k), empty(m,k), wait(m,k)

% FORMULATE PROBLEM

% All slots are integers
IntVars = ones(n,1);
x_L = zeros(n,1);
x_U = [ones(n1*n2,1);inf*ones(n3*(n1-1) + (n3-1)*n2,1)];

% Assignment constraint one job per rank.
A1 = zeros(n1,n);
for i=1:n1
   A1(i,(i-1)*n1+1:i*n1) = ones(1,n2);
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% Assignment constraint one rank per job.
A2 = zeros(n2,n);
for i=1:n2
   A2(i,i:n1:n2*n1-n1+i) = ones(1,n1);
end
b_L2 = ones(n2,1);
b_U2 = ones(n2,1);

% Relationship between the end of job rank k on machine m and start of job
% on machines m+1

A3 = zeros((n2-1)*(n3-1),n);
counter = 1;
for m=1:n3-1
   for k=1:n2-1
      A3(counter,[n1*n2+m+n3*(k-1), k*n1+1:(k+1)*n1, n1*n2+(n1-1)*n3+m+(n3-1)*k, ...
         n1*n2+(n1-1)*n3+m+(n3-1)*(k-1), (k-1)*n1+1:k*n1, n1*n2+m+1+n3*(k-1)]) = ...
         [1, proctimes(m,:), 1, -1, -proctimes(m+1,:), -1];
      counter = counter + 1;
   end
end
b_L3 = zeros((n2-1)*(n3-1),1);
b_U3 = zeros((n2-1)*(n3-1),1);

% Empty and Wait are zeros when starting, set in bounds
x_U(n1*n2+1:n3:n1*n2+(n1-1)*n3)             = zeros(n2-1,1);
x_U(n1*n2+(n1-1)*n3+1:n1*n2+(n1-1)*n3+n3-1) = zeros(n3-1,1);

% Add A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
c = zeros(n,1);
c(1:n1,1) = sum(proctimes(1:end-1,:),1)';
c(n1*n2+n3:n3:n1*n2+n3*(n1-1),1) = ones(n1-1,1); 

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Flow Shop Scheduling', [], [], IntVars);

% MODIFICATION LOG
%
% 051010 med   Created.