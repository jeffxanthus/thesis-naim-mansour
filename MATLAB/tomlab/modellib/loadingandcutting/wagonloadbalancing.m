% function Prob = wagonloadbalancing(minload, maxload, weights)
%
% Creates a TOMLAB MIP problem for wagon loading balancing (to be solved
% with infLinSolve)
%
% INPUT PARAMETERS
% minload      Minimum load for a wagon.
% maxload      Maximum load for the wagon.
% weights      Weight of the boxes.
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type LP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Prob = wagonloadbalancing(minload, maxload, weights)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(minload) | isempty(maxload) | isempty(weights)
   error('One of the inputs are empty');
end

m = length(minload);

if m~=length(maxload)
   error('Incorrect sizes for either minload or maxload');
end

minload = minload(:);
maxload = maxload(:);
weights = weights(:);

% FORMULATE PROBLEM
n1 = length(weights);
n = m*n1;

c = zeros(n,1); % Dummy objective
x_L = zeros(n,1);
x_U = ones(n,1);
IntVars = ones(n,1);

% Load constraint, exactly one box per wagon

A1 = repmat(speye(n1),1,m);
b_L1 = ones(n1,1);
b_U1 = b_L1;

% Wagon constraint

A2 = zeros(m,n);
for i=1:m
   A2(i,(i-1)*n1+1:i*n1) = weights';
end

b_L2 = minload;
b_U2 = maxload;

% Production and component constraint together

A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Wagon Load Balancing',...
   [], [], IntVars);

% Defining objective for minimax

Prob.QP.D = A2;

% MODIFICATION LOG
%
% 051007 med   Created.
% 060111 per   Changed line 49, ch 25 to 'm' instead of '3' - now more generic