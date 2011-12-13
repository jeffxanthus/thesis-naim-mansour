% function Prob = financinganearlyretirement(amounts, baseinterest,
% values, interest, duration)
%
% Creates a TOMLAB MIP problem for financing an early retirement scheme
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Prob = financinganearlyretirement(amounts, baseinterest, values, interest, duration)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(amounts) | isempty(baseinterest) | isempty(values) ...
      | isempty(interest) | isempty(duration)
   error('One of the inputs are empty');
end

n1  = length(values);    % Bonds to buy
n2  = length(amounts)-1; % Invest occasions 
n3  = 1;                 % Capital
n   = n1+n2+n3;
% FORMULATE PROBLEM
% Some variables are integer.
IntVars   = [ones(n1,1);zeros(n2+n3,1)];
x_L       = zeros(n,1);
x_U       = inf*ones(n,1);

% First year constr
A1   = [-values', -1, zeros(1,n2-1), ones(1,n3)];
b_L1 = amounts(1);
b_U1 = amounts(1);

% 2nd - 4th year constr.
A2   = zeros(3,n);
for i=1:3
   A2(i,[1:n1, n1+i:n1+i+1]) = [values'.*interest'/100, (1+baseinterest/100), -1];
end
b_L2 = amounts(2:4);
b_U2 = amounts(2:4);

% 5th - 6th year constr.
A3   = zeros(2,n);
for i=1:2
   idx1 = find(duration == i+3);
   idx2 = find(duration >= i+4);
   A3(i,[idx1', idx2', n1+i+3:n1+i+4]) = [values(idx1)'.*(1+interest(idx1)'/100), ...
      values(idx2)'.*interest(idx2)'/100, (1+baseinterest/100), -1];
end
b_L3 = amounts(5:6);
b_U3 = amounts(5:6);

% 7th year constr.
A4   = zeros(1,n);
idx = find(duration == 6);
A4(1,[idx, n1+6]) = [values(idx)'.*(1+interest(idx)'/100), (1+baseinterest/100)];
b_L4 = amounts(7);
b_U4 = amounts(7);

A   = [A1;A2;A3;A4];
b_L = [b_L1;b_L2;b_L3;b_L4];
b_U = [b_U1;b_U2;b_U3;b_U4];
c   = zeros(n,1);
c(end,1) = 1;
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Financing an Early Retirement Scheme', [], [], IntVars);

% MODIFICATION LOG
%
% 0151201 med   Created.