% function Prob = cuttingsheetmetal(demand, patterns)
%
% Creates a TOMLAB MIP problem for cutting smaller parts from big ones.
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = cuttingsheetmetal(demand, patterns)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(demand) | isempty(patterns)
   error('One of the inputs are empty');
end

demand   = demand(:);

if length(demand) ~= size(patterns,1)
   error('Size of pattern doesn''t match demand');
end

% FORMULATE PROBLEM
m = size(patterns,1); % cutting patterns
n = size(patterns,2);

c = ones(n,1);
x_L = zeros(n,1);
x_U = inf*ones(n,1);
IntVars = ones(n,1);

% Minimum demand must be met
A = patterns;
b_L = demand;
b_U = inf*ones(m,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Cutting Sheet Metal',...
   [],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.