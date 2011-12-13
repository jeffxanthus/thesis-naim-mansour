% function Prob = cuttingsteelbarsfordesklegs(demand, patterns, lengths)
%
% Creates a TOMLAB MIP problem for cutting steel bars for desk legs.
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = cuttingsteelbarsfordesklegs(demand, patterns, lengths)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(demand) | isempty(patterns) | isempty(lengths)
   error('One of the inputs are empty');
end

demand   = demand(:);
lengths  = lengths(:);

if length(demand) ~= size(patterns,1)
   error('Size of pattern doesn''t match demand');
end

% FORMULATE PROBLEM
m = size(patterns,1); % cutting patterns
n = size(patterns,2);

n2 = length(lengths);
n1 = n/length(lengths);

if abs(round(n1)-n1) > 1e-6
   error('pattern length have be dividable by lengths');
end

c = zeros(n,1);
for i=1:n2
   c((i-1)*n1+1:i*n1,1) = lengths(i);
end
x_L = zeros(n,1);
x_U = inf*ones(n,1);
IntVars = ones(n,1);

% Minimum demand must be met
A = patterns;
b_L = demand;
b_U = inf*ones(m,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Cutting Steel Bars',...
   [],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.