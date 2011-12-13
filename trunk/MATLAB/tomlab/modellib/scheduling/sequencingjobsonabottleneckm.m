% function Prob = sequencingjobsonabottleneckm(releasedate, duration,
% duedate, idx)
%
% Creates a TOMLAB MIP problem for sequencing jobs on a bottleneck machine
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Jan 2, 2005.   Last modified Jan 2, 2006.

function Prob = sequencingjobsonabottleneckm(releasedate, duration, duedate, idx)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(releasedate) | isempty(duration) | isempty(duedate) | isempty(idx)
   error('One of the inputs are empty');
end

if idx == 1

   n1 = length(releasedate); % Number of jobs     (j in JOBS)
   n2 = n1*n1;
   n  = n2+n1;

   % FORMULATE PROBLEM
   % All slots are integers
   IntVars = ones(n,1);
   x_L = zeros(n,1);
   x_U = [ones(n2,1);inf*ones(n1,1)];

   % Assignment constraint one job per rank.
   A1 = zeros(n1,n);
   for i=1:n1
      A1(i,(i-1)*n1+1:i*n1) = ones(1,n1);
   end
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);

   % Assignment constraint one rank per job.
   A2 = zeros(n1,n);
   for i=1:n1
      A2(i,i:n1:n2-n1+i) = ones(1,n1);
   end
   b_L2 = ones(n1,1);
   b_U2 = ones(n1,1);
   
   % Release constraints.
   A3 = zeros(n1,n);
   for i=1:n1
      A3(i,[n2+i, (i-1)*n1+1:i*n1]) = [1, -releasedate'];
   end
   b_L3 = zeros(n1,1);
   b_U3 = inf*ones(n1,1);
   
   % No simultaneous constraints
   A4 = zeros(n1-1,n);
   for i=1:n1-1
      A4(i,[n2+i+1, n2+i, (i-1)*n1+1:i*n1]) = [1, -1, -duration'];
   end
   b_L4 = zeros(n1-1,1);
   b_U4 = inf*ones(n1-1,1);
   
   % Add A, b_L, b_U
   A = [A1;A2;A3;A4];
   b_L = [b_L1;b_L2;b_L3;b_L4];
   b_U = [b_U1;b_U2;b_U3;b_U4];

   % Objective
   c = zeros(n,1);
   c(end,1) = 1;
   c(n2-n1+1:n2,1) = duration;
end

if idx == 2

   n1 = length(releasedate); % Number of jobs     (j in JOBS)
   n2 = n1*n1;
   n  = n2+n1+n1;

   % FORMULATE PROBLEM
   % All slots are integers
   IntVars = ones(n,1);
   x_L = zeros(n,1);
   x_U = [ones(n2,1);inf*ones(n1+n1,1)];

   % Assignment constraint one job per rank.
   A1 = zeros(n1,n);
   for i=1:n1
      A1(i,(i-1)*n1+1:i*n1) = ones(1,n1);
   end
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);

   % Assignment constraint one rank per job.
   A2 = zeros(n1,n);
   for i=1:n1
      A2(i,i:n1:n2-n1+i) = ones(1,n1);
   end
   b_L2 = ones(n1,1);
   b_U2 = ones(n1,1);
   
   % Release constraints.
   A3 = zeros(n1,n);
   for i=1:n1
      A3(i,[n2+i, (i-1)*n1+1:i*n1]) = [1, -releasedate'];
   end
   b_L3 = zeros(n1,1);
   b_U3 = inf*ones(n1,1);
   
   % No simultaneous constraints
   A4 = zeros(n1-1,n);
   for i=1:n1-1
      A4(i,[n2+i+1, n2+i, (i-1)*n1+1:i*n1]) = [1, -1, -duration'];
   end
   b_L4 = zeros(n1-1,1);
   b_U4 = inf*ones(n1-1,1);
   
   % Completion time equality
   A5 = zeros(n1,n);
   for i=1:n1
      A5(i,[n2+n1+i, n2+i, (i-1)*n1+1:i*n1]) = [1, -1, -duration'];
   end
   b_L5 = zeros(n1,1);
   b_U5 = zeros(n1,1);
   
   % Add A, b_L, b_U
   A = [A1;A2;A3;A4;A5];
   b_L = [b_L1;b_L2;b_L3;b_L4;b_L5];
   b_U = [b_U1;b_U2;b_U3;b_U4;b_U5];

   % Objective
   c = zeros(n,1);
   c(n2+n1+1:n,1) = ones(n1,1);
end

if idx == 3

   n1 = length(releasedate); % Number of jobs     (j in JOBS)
   n2 = n1*n1;
   n  = n2+n1+n1+n1;

   % FORMULATE PROBLEM
   % All slots are integers
   IntVars = ones(n,1);
   x_L = zeros(n,1);
   x_U = [ones(n2,1);inf*ones(n1+n1+n1,1)];

   % Assignment constraint one job per rank.
   A1 = zeros(n1,n);
   for i=1:n1
      A1(i,(i-1)*n1+1:i*n1) = ones(1,n1);
   end
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);

   % Assignment constraint one rank per job.
   A2 = zeros(n1,n);
   for i=1:n1
      A2(i,i:n1:n2-n1+i) = ones(1,n1);
   end
   b_L2 = ones(n1,1);
   b_U2 = ones(n1,1);
   
   % Release constraints.
   A3 = zeros(n1,n);
   for i=1:n1
      A3(i,[n2+i, (i-1)*n1+1:i*n1]) = [1, -releasedate'];
   end
   b_L3 = zeros(n1,1);
   b_U3 = inf*ones(n1,1);
   
   % No simultaneous constraints
   A4 = zeros(n1-1,n);
   for i=1:n1-1
      A4(i,[n2+i+1, n2+i, (i-1)*n1+1:i*n1]) = [1, -1, -duration'];
   end
   b_L4 = zeros(n1-1,1);
   b_U4 = inf*ones(n1-1,1);
   
   % Completion time equality
   A5 = zeros(n1,n);
   for i=1:n1
      A5(i,[n2+n1+i, n2+i, (i-1)*n1+1:i*n1]) = [1, -1, -duration'];
   end
   b_L5 = zeros(n1,1);
   b_U5 = zeros(n1,1);
   
   % Tardiness constraint
   A6 = zeros(n1,n);
   for i=1:n1
      A6(i,[n2+n1+n1+i, n2+n1+i, (i-1)*n1+1:i*n1]) = [1, -1, duedate'];
   end
   b_L6 = zeros(n1,1);
   b_U6 = inf*ones(n1,1);
   
   % Add A, b_L, b_U
   A = [A1;A2;A3;A4;A5;A6];
   b_L = [b_L1;b_L2;b_L3;b_L4;b_L5;b_L6];
   b_U = [b_U1;b_U2;b_U3;b_U4;b_U5;b_U6];

   % Objective
   c = zeros(n,1);
   c(n2+n1+n1+1:n,1) = ones(n1,1);
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Sequencing with Bottleneck', [], [], IntVars);

% MODIFICATION LOG
%
% 060102 med   Created.