% function Prob = productionofelectricity(demand, available, mincap, maxcap,...
%                 fixcost, runningcost, startcost);
%
% Creates a TOMLAB LP problem for production of electricity
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Prob = productionofelectricity(demand, available, mincap, maxcap,...
                 fixcost, runningcost, startcost, periodlengths)

if nargin < 8
   error('The function requires 8 inputs');
end

if isempty(demand) | isempty(available) | isempty(mincap) |...
      isempty(maxcap) | isempty(fixcost) | isempty(runningcost) |...
      isempty(startcost) | isempty(periodlengths)
   error('One of the inputs are empty');
end

demand       = demand(:);
available    = available(:);
mincap       = mincap(:);
maxcap       = maxcap(:);
fixcost      = fixcost(:);
runningcost  = runningcost(:);
startcost    = startcost(:);
periodlengths = periodlengths(:);

n1 = length(demand);
n2 = length(available);
n = n1*(3*n2);

% start_pt, work_pt, padd_pt

% FORMULATE PROBLEM

% All slots are integers
IntVars = ones(n,1);
IntVars(n-n1*n2+1:n) = 0;
x_L = zeros(n,1);
addcap = repmat((maxcap-mincap),n1,1);
x_U = [repmat(available,2*n1,1);inf*ones(length(addcap),1)];

% Relationship between work_pt and padd_pt
m1 = n1*n2;
A1 = zeros(m1,n);
for i=1:m1
   A1(i,2*m1+i) = 1;
   A1(i,m1+i)   = -addcap(i);
end
b_L1 = -inf*ones(m1,1);
b_U1 = zeros(m1,1);

% Demand in each period constraint
m2 = n1;
A2 = zeros(m2,n);
for i=1:m2
   A2(i,n1*n2+(i-1)*n2+1:n1*n2+i*n2) = mincap';
   A2(i,2*n1*n2+(i-1)*n2+1:2*n1*n2+i*n2) = ones(1,n2);
end
b_L2 = demand;
b_U2 = inf*ones(m2,1);

% 20% additional production possible
m3 = n1;
A3 = zeros(m3,n);
for i=1:m3
   A3(i,m1+(i-1)*n2+1:m1+i*n2) = maxcap';
end
b_L3 = 1.2*demand;
b_U3 = inf*ones(m3,1);

% Relationship between start_pt and work_pt

m4 = n1*n2;
A4 = zeros(m4,n);
for i=1:n2
   A4(i,i) = 1;
   A4(i,n1*n2+i) = -1;
   A4(i,2*n1*n2-n2+i) = 1;
end
for i=n2+1:m4
   A4(i,i) = 1;
   A4(i,n1*n2+i) = -1;
   A4(i,n1*n2+i-n2) = 1;
end

b_L4 = zeros(m4,1);
b_U4 = inf*ones(m4,1);

% Add A, b_L, b_U

A = [A1;A2;A3;A4];
b_L = [b_L1;b_L2;b_L3;b_L4];
b_U = [b_U1;b_U2;b_U3;b_U4];

% Objective

c = zeros(n,1);
c(1:n1*n2,1) = repmat(startcost,n1,1);
vec1 = fixcost*periodlengths';
vec1 = vec1(:);
c(n1*n2+1:2*n1*n2,1) = vec1;
vec2 = runningcost*periodlengths';
vec2 = vec2(:);
c(2*n1*n2+1:n,1) = vec2;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Production of Electricity', [], [], IntVars);

% MODIFICATION LOG
%
% 051010 med   Created.