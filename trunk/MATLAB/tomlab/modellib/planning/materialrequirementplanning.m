% function Prob = materialrequirementplanning(demand, compprices, assemmat...
%   , subcontr, assembly, finalassembly, capacity, finalcapacity)
%
% Creates a TOMLAB MIP problem for material requirement planning
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 18, 2005.   Last modified Oct 18, 2005.

function Prob = materialrequirementplanning(demand, compprices, assemmat...
   , subcontr, assembly, finalassembly, capacity, finalcapacity)

if nargin < 8
   error('The function requires 8 inputs');
end

if isempty(demand) | isempty(compprices) | isempty(assemmat)...
      | isempty(subcontr) | isempty(assembly) ...
      | isempty(finalassembly) | isempty(capacity) | isempty(finalcapacity)
   error('One of the inputs are empty');
end

n1 = length(compprices); % Buy preprod (12), Buy subcontr (3), assemble (3), finalass (2)
n2 = length(subcontr);
n3 = length(assembly);
n4 = length(finalassembly);
n = n1+n2+n3+n4;

% FORMULATE PROBLEM

% All slots are integers
IntVars    = zeros(n,1);
x_L        = [zeros(n1+n2+n3,1);demand];
x_U        = [inf*ones(n1+n2,1);capacity;finalcapacity];

% For assembled products, some has to be 2 times the other
A   = assemmat;
b_L = zeros(size(A,1),1);
b_U = inf*ones(size(A,1),1);

% Objective
c = [compprices;subcontr;assembly;finalassembly];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Material Requirement Planning', [], [], IntVars);

% MODIFICATION LOG
%
% 051018 med   Created.