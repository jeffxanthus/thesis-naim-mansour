% function Prob = cctvsurveillance(arcs_in, arcs_out)
%
% Creates a TOMLAB MILP problem for cctv surveillance
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = cctvsurveillance(arcs_in, arcs_out)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(arcs_in) | isempty(arcs_out)
   error('One of the inputs are empty');
end

n    = length(arcs_in);  %arcs

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% All streets need to be covered
b_L  = ones(n,1);
b_U  = inf*ones(n,1);
A    = zeros(n,n);
for i=1:n
   A(i, [arcs_in(i), arcs_out(i)]) = [1 1];
end

c = ones(n,1); 

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'CCTV Surveillance', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.