% function Prob = canesugarproduction(proclines, proctime, loss, lifespan)
%
% Creates a TOMLAB LP problem for cane sugar production
%
% INPUT PARAMETERS
% proclines     Number of production lines
% proctime      Processing time (generic unit 2).
% loss          Loss (in generic unit 1 / generic unit 2).
% lifespan      Life span of . (in generic unit 1).
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Prob = canesugarproduction(proclines, proctime, loss, lifespan)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(proctime) | isempty(loss) | isempty(lifespan) | isempty(proclines) | isempty(proctime)
    error('One of the inputs are empty');
end

loss     = loss(:);   
lifespan = lifespan(:);

n = length(loss);

if n~=length(lifespan)
   error('Incorrect sizes for either loss or lifespan');
end

% FORMULATE PROBLEM

n1 = length(loss);
n2 = ceil(n1/proclines);
n   = n1*n2;

% All slots are integers
IntVars = ones(n,1);
x_L = zeros(n,1);
x_U = ones(n,1);

% Production constraint, all lots processed
A1 = zeros(n1,n);
for i=1:n1
   A1(i,i:n1:n-n1+i) = 1;
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% Production constraint, No more than proclines concurrently running

A2 = zeros(n2,n);
for i=1:n2
   A2(i,(i-1)*n1+1:i*n1) = 1;
end
b_L2 = -inf*ones(n2,1);
b_U2 = ones(n2,1)*proclines;

% Maximum slot number for a lot
vec = repmat([1:n2]',1,n1);
vec = vec(:);
A3 = zeros(n1,n);
for i=1:n1
   A3(i,i:n1:n-n1+i) = vec((i-1)*n2+1:i*n2)';
end
b_L3 = -inf*ones(n1,1);
b_U3 = lifespan/proctime;

A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

loss = repmat(loss,n2,1);
vec = repmat([1:n2]',1,n1);
vec = vec';
vec = vec(:);
c = vec.*loss*proctime;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Cane Sugar Production', [], [], IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.