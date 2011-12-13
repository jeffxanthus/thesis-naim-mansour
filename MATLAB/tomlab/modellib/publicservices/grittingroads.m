% function Prob = grittingroads()
%
% Creates a TOMLAB MILP problem for gritting roads
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 6, 2005.   Last modified Dec 6, 2005.

function Prob = grittingroads(in, out, lengths)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(in) | isempty(out) | isempty(lengths)
   error('One of the inputs are empty');
end

n    = length(lengths);     %Number of arcs
n1   = length(unique([in;out]));

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = ones(n,1);
x_U       = inf*ones(n,1);

% Districts chosen
b_L  = zeros(n1,1);
b_U  = zeros(n1,1);
A    = zeros(n1,n);
for i=1:n1
   idxin  = find(i == in);
   idxout = find(i == out);
   A(i,[idxout;idxin]) = [ones(1,length(idxout)), -ones(1,length(idxin))];
end

c = lengths; 

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Gritting Roads', [], [], IntVars);

% MODIFICATION LOG
%
% 051206 med   Created.