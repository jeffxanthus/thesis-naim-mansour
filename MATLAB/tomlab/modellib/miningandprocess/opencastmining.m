% function Prob = opencastmining(values, depends)
%
% Creates a TOMLAB LP problem for open cast mining
%
% INPUT PARAMETERS
% values        Profit for each extraction (generic unit 1).
% depends       Dependecies matrix (generic unit 2).
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = opencastmining(values, depends)

if nargin < 2
   error('The function requires 2 inputs');
end

if isempty(values) | isempty(depends)
   error('One of the inputs are empty');
end

values   = values(:);
n = length(values);

if n~=size(depends,2)
   error('Incorrect sizes for either values or depends');
end

% FORMULATE PROBLEM
m = size(depends,1);

% All slots are integers
IntVars = ones(n,1);
x_L = zeros(n,1);
x_U = ones(n,1);

% Production constraint, all lots processed
A = depends;
b_L = -inf*ones(m,1);
b_U = zeros(m,1);

c = -values;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Open Cast Mining', [], [], IntVars);

% MODIFICATION LOG
%
% 051010 med   Created.