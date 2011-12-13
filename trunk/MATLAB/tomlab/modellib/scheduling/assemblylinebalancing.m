% function Prob = assemblylinebalancing(stations, duration, dependsvec1,
% dependsvec2)
%
% Creates a TOMLAB MIP problem for assembly line balancing
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 13, 2005.   Last modified Jan 2, 2006.

function Prob = assemblylinebalancing(stations, duration, dependsvec1, dependsvec2)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(stations) | isempty(duration) | isempty(dependsvec1) | isempty(dependsvec2)
   error('One of the inputs are empty');
end

n1 = length(duration);
n2 = stations;
n  = n1*n2 + 1; % Final variables is cycle time

duration = duration(:);
dependsvec1 = dependsvec1(:);
dependsvec2 = dependsvec2(:);

m = length(dependsvec1);

% FORMULATE PROBLEM

% All slots are integers
IntVars    = ones(n-1,1);
x_L        = zeros(n,1);
x_U        = ones(n,1);
x_U(end,1) = inf;

% Every process has to go to one machine
A1 = zeros(n1,n);
for i=1:n1
   A1(i,i:n1:n-1-n1+i) = 1;
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% Precedence constraint
A1 = zeros(n1,n);
for i=1:n1
   A1(i,i:n1:n-1-n1+i) = 1;
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

A2 = zeros(m,n);
for i=1:m
   idx1 = [dependsvec1(i):n1:n1*n2-n1+dependsvec1(i)]';
   idx2 = [dependsvec2(i):n1:n1*n2-n1+dependsvec2(i)]';
   A2(i,idx1) = -[1:n2];
   A2(i,idx2) = [1:n2];
end

b_L2 = -inf*ones(m,1);
b_U2 = zeros(m,1);

% Cycle constraint
A3 = zeros(n2,n);
for i=1:n2
   A3(i,(i-1)*n1+1:i*n1) = duration';
   A3(i,end) = -1; % Cycle
end

b_L3 = -inf*ones(n2,1);
b_U3 = zeros(n2,1);

% Add A, b_L, b_U
A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
c = zeros(n,1);
c(end,1) = 1;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Assembly Line Balancing', [], [], IntVars);

% MODIFICATION LOG
%
% 051017 med   Created
% 060102 med   Model verified