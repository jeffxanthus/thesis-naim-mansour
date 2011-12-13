% function Prob = composingflightcrews(scores, arcs_in, arcs_out,
% arcs_score)
%
% Creates a TOMLAB MIP problem for composing flight crews
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 21, 2005.   Last modified Oct 21, 2005.

function Prob = composingflightcrews(scores, arcs_in, arcs_out, arcs_score)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(scores) | isempty(arcs_in) | isempty(arcs_out) | isempty(arcs_score)
   error('One of the inputs are empty');
end

n = length(arcs_in);         % psbl
n1 = size(scores,2);
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = ones(n,1);

% orig constr.
A = zeros(n1,n);
for i=1:n1
   idx_in  = find(arcs_in  == i);
   idx_out = find(arcs_out == i);
   idx = unique([idx_in;idx_out]);
   A(i,idx) = ones(1,length(idx));
end
b_L = -inf*ones(n1,1);
b_U = ones(n1,1);

% Objective
c   = -arcs_score;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Composing Flight Crews', [], [], IntVars);

% MODIFICATION LOG
%
% 051021 med   Created.