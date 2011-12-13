% function Prob = waterconveyance(arcs_in, arcs_out, capacity, source,
% sink)
%
% Creates a TOMLAB MILP problem for water conveyance
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = waterconveyance(arcs_in, arcs_out, capacity, source, sink)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(arcs_in) | isempty(arcs_out) | isempty(capacity) | isempty(source)...
      | isempty(sink)
   error('One of the inputs are empty');
end

n    = length(arcs_in);  %arcs

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = capacity;

% Kirchhoffs's law
n1   = length(unique([arcs_in;arcs_out])) - 2;
b_L  = zeros(n1,1);
b_U  = zeros(n1,1);
A    = zeros(n1,n);
for i=1:n1+2
   if ~(i ==  source | i == sink)
      idx1 = find(arcs_in == i);
      idx2 = find(arcs_out == i);
      A(i, [idx1;idx2]) = [ones(1,length(idx1)), -ones(1,length(idx2))];
   end
end

idx = find(arcs_out == sink);
c   = zeros(n,1);
c(idx,1) = -ones(length(idx),1); 

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Water Conveyance', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.