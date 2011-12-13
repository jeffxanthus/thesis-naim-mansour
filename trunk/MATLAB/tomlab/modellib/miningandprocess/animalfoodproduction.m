% function Prob = animalfoodproduction(compsize, mincomp, maxcomp, rawcompmat,
% rawavail, rawcost, prodcostsraw, prodcostsprod)
%
% Creates a TOMLAB LP problem for animal food production
%
% INPUT PARAMETERS
% compsize      Size of component order (in generic unit 2).
% mincomp       Lower limit on specification. Vector. (in generic unit 1).
% maxcomp       Upper limit on specification. Vector. (in generic unit 1).
%
% rawcompmat    Component matrix for raw materials (in generic unit 1).
% rawavail      Availability (in general unit 2).
% rawcost       Cost (in generic currency / general unit 2)
% prodcostsraw  Raw material processing costs.
% prodcostsprod Prodcu material processing costs.
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type LP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Prob = animalfoodproduction(compsize, mincomp, maxcomp, rawcompmat,...
   rawavail, rawcost, prodcostsraw, prodcostsprod)

if nargin < 7
   error('The function requires 7 inputs');
end

if isempty(mincomp) | isempty(maxcomp) | isempty(rawcompmat) |...
       isempty(rawavail) | isempty(rawcost) | isempty(prodcostsraw) | isempty(prodcostsprod)
    error('One of the inputs are empty');
end

mincomp = mincomp(:); % lower bounds on components
maxcomp = maxcomp(:); % upper bounds on components

m = length(mincomp);

if m~=length(maxcomp) | m~=size(rawcompmat,2)
   error('Incorrect sizes for either maxcomp, rawcompmat or prodcosts (columns should equal rows in maxcomp');
end

rawavail = rawavail(:);
rawcost  = rawcost(:);

n = size(rawcompmat,1);

if n~=length(rawavail) | n~=length(rawcost)
   error('Incorrect sizes for either rawavail or rawcost');
end

compsize = compsize(:);

% FORMULATE PROBLEM

% n is raw * prods + prods
n1  = n;
n2  = size(prodcostsprod,2);
n   = n*n2+n2;

% Cost of raw material
c   = [repmat(rawcost,n2,1);zeros(n2,1)];

% Cost of production
c  = c + sum([repmat(prodcostsraw,1,n2),prodcostsprod],1)';
x_L = [zeros(n-n2,1);compsize];
x_U = [repmat(rawavail,n2,1);compsize];

% Production constraint
A1 = zeros(n2,n);
for i=1:n2
   A1(i,(i-1)*n1+1:i*n1) = ones(1,n1);
end
b_L1 = compsize;
b_U1 = compsize;

% Component constraint
A2   = zeros(n2*size(rawcompmat',1),n);
b_L2 = repmat(mincomp,n2,1);
b_U2 = repmat(maxcomp,n2,1);

for i=1:n2
   A2((i-1)*size(rawcompmat',1)+1:i*size(rawcompmat',1),(i-1)*n1+1:i*n1) = rawcompmat';
   b_L2((i-1)*length(mincomp)+1:i*length(mincomp)) = b_L2((i-1)*length(mincomp)+1:i*length(mincomp))*compsize(i);
   b_U2((i-1)*length(mincomp)+1:i*length(mincomp)) = b_U2((i-1)*length(mincomp)+1:i*length(mincomp))*compsize(i);
end

% Availability constraint
A3   = [eye(n1),eye(n1),zeros(n1,n2)];
b_L3 = -inf*ones(n1,1);
b_U3 = rawavail;

% Production and component constraint together

A = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, [], 'Animal Food Production');

% MODIFICATION LOG
%
% 051007 med   Created.