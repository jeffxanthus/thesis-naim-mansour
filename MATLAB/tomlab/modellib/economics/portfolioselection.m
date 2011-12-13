% function Prob = portfolioselection(budget, mininvest, maxinvest, catinvest1min, idx1cat,...
%   catinvest2max, idx2cat, returns)
%
% Creates a TOMLAB MIP problem for portfolio selection
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = portfolioselection(budget, mininvest, maxinvest, catinvest1min, idx1cat,...
   catinvest2max, idx2cat, returns)

if nargin < 8
   error('The function requires 8 inputs');
end

if isempty(budget) | isempty(mininvest) | isempty(maxinvest) | isempty(catinvest1min) ...
      | isempty(idx1cat) | isempty(catinvest2max) | isempty(idx2cat) | isempty(returns)
   error('One of the inputs are empty');
end

n  = length(returns);

% FORMULATE PROBLEM
% No variables are binary.
IntVars   = zeros(n,1);
x_L       = mininvest*ones(n,1);
x_U       = maxinvest*ones(n,1);

% Budget constraints
A1   = idx2cat';
b_L1 = -inf;
b_U1 = budget*catinvest2max;

A2   = idx1cat';
b_L2 = budget*catinvest1min;
b_U2 = inf;

A3   = ones(1,n);
b_L3 = budget;
b_U3 = budget;

A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

c   = -returns/100;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Portfolio Selection', [], [], IntVars);

% MODIFICATION LOG
%
% 051201 med   Created.