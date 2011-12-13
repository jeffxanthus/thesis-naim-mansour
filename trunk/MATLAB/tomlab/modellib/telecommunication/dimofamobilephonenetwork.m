% function Prob = dimofamobilephonenetwork(hub_mat, traffic,
% connections, capacity)
%
% Creates a TOMLAB MIP problem for dimensioning of a mobile phone network
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 8, 2005.   Last modified Nov 8, 2005.

function Prob = dimofamobilephonenetwork(hub_mat, traffic,...
   connections, capacity)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(hub_mat) | isempty(traffic) | isempty(connections) | isempty(capacity)
   error('One of the inputs are empty');
end

n1 = size(hub_mat,2);   % Cells
n2 = size(hub_mat,1);   % Nodes
n  = n1*n2;             % Cells (node1), Cells (node2)..., last node MTSO

% FORMULATE PROBLEM
% All variables are binary.
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Cells connected to minimun nodes
A1   = zeros(n1,n);
b_L1 = connections;
b_U1 = connections;
for i=1:n1
   A1(i,i:n1:n-n1+i) = ones(1,n2);
end

% Limits of the ring
A2   = zeros(1,n);
b_L2 = -inf*ones(1,1);
b_U2 = 2*capacity;
for i=1:n1 % For all cells
   for j=1:n2-1
      A2(1,(j-1)*n1+i) = traffic(i)/connections(i);
   end
end

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];
hub_mat = hub_mat';
c   = hub_mat(:);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Dimensioning of a Mobile Phone Network', [], [], IntVars);

% MODIFICATION LOG
%
% 051108 med   Created.