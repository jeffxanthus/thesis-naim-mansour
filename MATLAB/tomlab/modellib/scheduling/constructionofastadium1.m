% function Prob = constructionofastadium1(taskduration, taskprecedence);
%
% Creates a TOMLAB LP problem for production of a stadium
%
% INPUT PARAMETERS
% taskduration    Duration of tasks (in generic unit 1).
% taskprecedence  Matrix saying which tasks need to be completed before
%                 others (no unit).
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 10, 2005.   Last modified Jan 2, 2006.

function Prob = constructionofastadium1(taskduration, taskprecedence)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(taskduration) | isempty(taskprecedence)
    error('One of the inputs are empty');
end

taskduration = taskduration(:);

n = length(taskduration)+1; % Maximum final time slot

% FORMULATE PROBLEM

c   = [zeros(n-1,1);1];
x_L = zeros(n,1);
x_U = inf*ones(n,1);
IntVars = ones(n,1);

% Start constraint
counter = 1;
% Extra one for final time.
mLin = nnz(taskprecedence)+1; 
A1 = zeros(mLin,n);
b_L1 = -inf*ones(mLin,1);
b_U1  = zeros(mLin,1);

for i=1:n-1
   idx = find(taskprecedence(i,:) ~= 0);
   for j=1:length(idx)
      A1(counter,i) = -1;
      A1(counter,idx(j)) = 1;
      b_U1(counter,1) = -taskduration(i);
      counter = counter + 1;
   end
end
A1(end,end-1) =  1;
A1(end,end)   = -1;
b_L1(end,1)   = -inf;
b_U1(end,1)   = 0;

idx = find(sum(taskprecedence,2) == 0);
numidx = length(idx);
A2 = zeros(numidx,n);
b_L2 = -inf*ones(numidx,1);
b_U2 = zeros(numidx,1);
for i=1:length(idx)
   A2(i,i) = -1;
   b_U2(i) = -taskduration(idx);
end

A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Construction of a Stadium 1',[],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.