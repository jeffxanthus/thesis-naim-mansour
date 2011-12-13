% function Prob = familybudget(costs, costfreq, income, incomefreq, minhobby)
%
% Creates a TOMLAB MIP problem for family budget
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = familybudget(costs, costfreq, income, incomefreq, minhobby)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(costs) | isempty(costfreq) | isempty(income) | isempty(incomefreq)...
      | isempty(minhobby)
   error('One of the inputs are empty');
end

n1  = 12;                % Hobby, months
n2  = n1;                % Savings, months
n   = n1+n2;

% FORMULATE PROBLEM
% No variables are integer.
IntVars   = zeros(n,1);
x_L       = [minhobby*ones(n1,1);zeros(n2,1)];
x_U       = inf*ones(n,1);

% Months constraints
b_L = -inf*ones(n1,1);
b_U = zeros(n1,1);
value = 0;
incomefreqcount = incomefreq;
costfreqcount = costfreq;
for i=1:n1
   for j=1:length(incomefreq)
      if incomefreq(j) == i
         value = value + income(j);
         incomefreq(j) = incomefreq(j) + incomefreqcount(j);
      end
   end
   for j=1:length(costfreq)
      if costfreq(j) == i
         value = value - costs(j);
         costfreq(j) = costfreq(j) + costfreqcount(j);
      end
   end
   b_U(i,1) = value;
   value = 0;
end

A   = zeros(n1,n);
A(1,[1,1+n1]) = [1 1];
for i=2:n1
   A(i,[i,i-1+n1,i+n1]) = [1 -1 1];
end

c   = -[ones(n1,1);zeros(n2,1)];
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Family Budget', [], [], IntVars);

% MODIFICATION LOG
%
% 0151201 med   Created.