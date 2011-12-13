% function Prob = bargeloading(capacity, units, sizes, unitprice, cost)
%
% Creates a TOMLAB MIP problem for barge loading.
%
% INPUT PARAMETERS
% capacity     Capacity of the barge
% units        Number of units available from each location
% sizes        Size of units at locations (generic unit 1)
% unitprice    Price charged for each location (generic currency)
% cost         Cost to ship (generic currency / generic unit 1)
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Prob = bargeloading(capacity, units, sizes, unitprice, cost)

if nargin < 5
   error('The function requires 5 inputs');
end

if isempty(capacity) | isempty(units) | isempty(sizes) | isempty(unitprice) | isempty(cost)
   error('One of the inputs are empty');
end

m = length(units);

if m~=length(sizes) | m~=length(unitprice) | m~=length(cost)
   error('Incorrect sizes for one of the inputs');
end

units = units(:);
sizes = sizes(:);
unitprice = unitprice(:);
cost = cost(:);

% FORMULATE PROBLEM
n = length(units);

c = -(unitprice - sizes.*cost); % - since maximising
x_L = zeros(n,1);
x_U = units;
IntVars = ones(n,1);

% Load constraint, cannot exceed capacity

A    = sizes';
b_L  = capacity;
b_U  = b_L;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Barge Loading',...
   [],[],IntVars);

% MODIFICATION LOG
%
% 051007 med   Created.