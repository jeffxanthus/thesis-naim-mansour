% function Prob = choiceofexpansionprojects(benefit, budget, costmat)
%
% Creates a TOMLAB MIP problem for choice of expansion projects
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = choiceofexpansionprojects(benefit, budget, costmat)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(benefit) | isempty(budget) | isempty(costmat)
   error('One of the inputs are empty');
end

n1  = size(costmat,2);   %years 
n   = length(benefit);   %projects

% FORMULATE PROBLEM
% All variables are integer.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Cost constraints
b_L = -inf*ones(n1,1);
b_U = budget;
A   = costmat';

c   = -benefit;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Choice of Expansion Projects', [], [], IntVars);

% MODIFICATION LOG
%
% 0151201 med   Created.