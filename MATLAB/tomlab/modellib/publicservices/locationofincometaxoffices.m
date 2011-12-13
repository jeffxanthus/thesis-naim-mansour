% function Prob = locationofincometaxoffices(population, numloc, in, out,
% lengths)
%
% Creates a TOMLAB MILP problem for location of income tax offices
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 6, 2005.   Last modified Dec 6, 2005.

function Prob = locationofincometaxoffices(population, numloc, in, out, lengths)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(population) | isempty(numloc) | isempty(in) | isempty(out) | isempty(lengths)
   error('One of the inputs are empty');
end

n1   = length(unique([in;out])); %Number of cities
n2   = length(in);
% Calculate distance matrix

d = inf*ones(n1,n1);
for i=1:n1
   d(i,i) = 0;
end
for i=1:n2
   d(in(i), out(i)) = lengths(i);
   d(out(i), in(i)) = lengths(i);
end
for i=1:n1 %b
   for j=1:n1 %c
      for k=1:n1 %d
         if j<k
            if d(j,k) > d(j,i)+d(i,k);
               d(j,k) = d(j,i)+d(i,k);
               d(k,j) = d(j,i)+d(i,k);
            end
         end
      end
   end
end

n = n1+n1^2;

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Building constraint
b_L1  = numloc;
b_U1  = numloc;
A1    = [ones(1,n1), zeros(1,n1^2)];

% Dependencies constraint
b_L2  = ones(n1,1);
b_U2  = ones(n1,1);
A2    = zeros(n1,n);
for i=1:n1
   A2(i,n1+(i-1)*n1+1:n1+i*n1) = ones(1,n1);
end

% Reality constraint
b_L3  = -inf*ones(n-n1,1);
b_U3  = zeros(n-n1,1);
A3    = zeros(n-n1,n);
for i=1:n1
   for j=1:n1
      A3((j-1)*n1+i,[i,n1+(j-1)*n1+i]) = [-1 1];
   end
end

% Merge constraints
A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

c = zeros(n,1);
for i=1:n1
   for j=1:n1
      c(n1+(i-1)*n1+j) = population(i)*d(i,j);
   end
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Location of Income Tax Offices', [], [], IntVars);

% MODIFICATION LOG
%
% 051206 med   Created.