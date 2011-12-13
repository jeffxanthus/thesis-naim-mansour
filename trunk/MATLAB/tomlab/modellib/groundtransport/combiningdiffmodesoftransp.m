% function Prob = combiningdiffmodesoftransp(transpcost,
% changecost, demand);
%
% Creates a TOMLAB MIP problem for combining different modes of transport
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 20, 2005.   Last modified Oct 20, 2005.

function Prob = combiningdiffmodesoftransp(transpcost, changecost, demand)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(transpcost) | isempty(changecost) | isempty(demand)
   error('One of the inputs are empty');
end

n1 = size(transpcost,2);   % legs
n2 = size(transpcost,1);   % mode
n3 = size(transpcost,2)-1; % change

n  = n1*n2 + n2*n2*n3;   % leg (mode 1), leg (mode 2)...
                         % change (mode-1*mode-2) change(mode-1*mode-2)...
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = ones(n,1);
x_L        = zeros(n,1);
x_U        = ones(n,1);

% Single mode of transport
A1 = zeros(n1,n);
for i=1:n1
   A1(i,(i-1)*n2+1:i*n2) = ones(1,n2);
end
b_L1 = ones(n1,1);
b_U1 = ones(n1,1);

% Single change of mode
A2 = zeros(n3,n);
for i=1:n3
   A2(i,n1*n2+(i-1)*n2*n2+1:n1*n2+i*n2*n2) = ones(1,n2*n2);
end
b_L2 = ones(n3,1);
b_U2 = ones(n3,1);

% Change to mode relationship
A3 = zeros(n2*n2*n3,n);
counter = 1;
for i=1:n2
   for j=1:n2
      for k=1:n3
         A3(counter,[(i-1)*n2+j, (i-1)*n2+n2+k ,n1*n2+counter ]) = [-1 -1 2];
         counter = counter + 1;
      end
   end
end
b_L3 = -inf*ones(n2*n2*n3,1);
b_U3 = zeros(n2*n2*n3,1);

% Merge A, b_L, b_U
A   = [A1;A2;A3];
b_L = [b_L1;b_L2;b_L3];
b_U = [b_U1;b_U2;b_U3];

% Objective
changecost = changecost';
c   = [transpcost(:);repmat(changecost(:),n3,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Combining Different Modes of Transport', [], [], IntVars);

% MODIFICATION LOG
%
% 051020 med   Created.