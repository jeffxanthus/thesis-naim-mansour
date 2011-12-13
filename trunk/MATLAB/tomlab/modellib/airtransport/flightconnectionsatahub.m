% function Prob = flightconnectionsatahub(origdest)
%
% Creates a TOMLAB MIP problem for flight connections at a hub
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 21, 2005.   Last modified Oct 21, 2005.

function Prob = flightconnectionsatahub(origdest, idximp)

if nargin < 1
   error('The function requires 1 inputs');
end

if isempty(origdest)
   error('One of the inputs are empty');
end

n1 = size(origdest,1);         % orig
n2 = size(origdest,2);         % dest

n  = n1*n2;              % dest, dest, dest

% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = ones(n,1);
x_U(idximp) = zeros(length(idximp),1);

% orig constr.
A1 = zeros(n1,n);
for i=1:n2
   A1(i,(i-1)*n2+1:i*n2) = ones(1,n1);
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% dest constr.
A2 = zeros(n2,n);
for i=1:n1
   A2(i,i:n1:n-n1+i) = ones(1,n2);
end
b_L2 = ones(n2,1);
b_U2 = ones(n2,1);

% Merge A, b_L, b_U

A   = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Objective
c   = -origdest(:);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Flight Connections at a Hub', [], [], IntVars);

% MODIFICATION LOG
%
% 051021 med   Created.