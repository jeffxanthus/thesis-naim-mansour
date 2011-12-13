% function Prob = heatingoildelivery(distance, demand, capacity);
%
% Creates a TOMLAB MIP problem for heating oil delivery
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 20, 2005.   Last modified Oct 20, 2005.

function Prob = heatingoildelivery(distance, demand, capacity)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(distance) | isempty(demand) | isempty(capacity)
   error('One of the inputs are empty');
end

n1 = size(distance,1);   % SITES
n2 = n1-1;               % CUSTOMERS
n3 = n1*n1;
n  = n3 + n2;            % sites (cust 1), sites (cust2)...
                         % quant (for all clients).
% FORMULATE PROBLEM

% All variables are integers.
IntVars    = [ones(n3,1);zeros(n2,1)];
x_L        = [zeros(n3,1);demand];
x_U        = [ones(n3,1);capacity*ones(n2,1)];
x_U(1:n1+1:n3,1) = zeros(n1,1);

% Customer constraint.
A1 = zeros(n2,n);
for j=2:n1
   for i=1:n1
      if j~=i
         A1(j-1,(j-1)*n1+i) = 1;
      end
   end
end
b_L1 = ones(n2,1);
b_U1 = ones(n2,1);

% Customer constraint.
A2 = zeros(n2,n);
for i=2:n1
   for j=1:n1
      if i~=j
         A2(i-1,i+(j-1)*n1) = 1;
      end
   end
end
b_L2 = ones(n2,1);
b_U2 = ones(n2,1);

% quantity, demand, capacity constr.
A3 = zeros(n2,n);
for i=2:n1
   A3(i-1,[n3+i-1, i]) = [1 -demand(i-1)+capacity];
end
b_L3 = -inf*ones(n2,1);
b_U3 = capacity*ones(n2,1);

% quantity, demand, capacity constr.
A4 = zeros(n2*n2-n2,n);
b_L4 = -inf*ones(n2*n2-n2,1);
b_U4 = zeros(n2*n2-n2,1);
counter = 1;
for i=2:n1
   for j=2:n1
      if i~=j
         A4(counter,[n3+i-1,n3+j-1,(i-1)*n1+j,i+(j-1)*n1]) = [1 -1 capacity capacity-demand(j-1)-demand(i-1)];
         b_U4(counter,1) = capacity-demand(j-1);
         counter = counter + 1;
      end
   end
end

% Merge A, b_L, b_U

A   = [A1;A2;A3;A4];
b_L = [b_L1;b_L2;b_L3;b_L4];
b_U = [b_U1;b_U2;b_U3;b_U4];

% Objective
c   = [distance(:);zeros(n2,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Heating Oil Delivery', [], [], IntVars);

% MODIFICATION LOG
%
% 051020 med   Created.