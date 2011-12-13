% function Prob = constructionofastadium2(maxreduc, costperweek, bonus,
% idx, taskduration, Result);
%
% Creates a TOMLAB MILP problem for project crashing
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 10, 2005.   Last modified Jan 2, 2006.

function Prob = constructionofastadium2(maxreduc, costperweek, bonus, idx, taskduration, Result)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(maxreduc) | isempty(costperweek) | isempty(bonus) | isempty(Result)...
      | isempty(taskduration) | isempty(idx)
   error('One of the inputs are empty');
end

maxreduc = maxreduc(:);
costperweek = costperweek(:);
bonus = bonus(:);

n1 = length(maxreduc);
n = n1*2+2; % Maximum final time slot

obj1 = Result.f_k;

% FORMULATE PROBLEM

% Save is first variable, start second
x_L = zeros(n,1);
x_U = [maxreduc;inf;Result.Prob.x_U];
IntVars = ones(n,1);

% Start-save constraint
A1 = zeros(1,n);
A1(1,n1+1) = 1;
A1(1,n)    = 1;
b_L1 = obj1;
b_U1 = obj1;

counter = 1;
A2   = [zeros(sum(idx)+2,n1+1),Result.Prob.A];
b_U2 = zeros(size(A2,1),1);
% Arc constraint
for i=1:length(idx)
   for j=1:idx(i)
      A2(counter,i) = -1;
      b_U2(counter,1) = -taskduration(i);
      counter = counter + 1;
   end
   if idx(i) == 0
      % Set bounds instead, no predecessors
      x_L(n1+1+i,1) = taskduration(i);
      x_U(n1+1+i,1) = taskduration(i);      
   end
end

b_L2 = -inf*ones(size(A2,1),1);

A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

c = [costperweek;-bonus;zeros(n1+1,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Construction of a Stadium 2',[],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.
% 060102 med   Fine-tuned problem