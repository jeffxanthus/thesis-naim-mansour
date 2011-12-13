% function Prob = planthepersonnelataconstr(transferin, transferout,...
%   staffingdevcost, overtimemax, maxtransferin, maxtransferout, ...
%   startstaff, endstaff, demands)
%
% Creates a TOMLAB MILP problem for planning the personnel at a construction site
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Prob = planthepersonnelataconstr(transferin, transferout,...
   staffingdevcost, overtimemax, maxtransferin, maxtransferout, ...
   startstaff, endstaff, demands)

if nargin < 9
   error('The function requires 9 inputs');
end

if isempty(transferin) | isempty(transferout) | isempty(staffingdevcost) | isempty(overtimemax)...
      | isempty(maxtransferin) | isempty(maxtransferout) | isempty(startstaff)...
      | isempty(endstaff) | isempty(demands)
   error('One of the inputs are empty');
end

n1    = length(demands);  %months
n     = 5*n1;             %onsite, arrive, leave, over, under

% FORMULATE PROBLEM
% No variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = inf*ones(n,1);

% Start constraint
b_L1 = startstaff;
b_U1 = startstaff;
A1   = zeros(1,n);
A1(1,[1, n1+1]) = [1 -1];

% Final constraint
b_L2 = endstaff;
b_U2 = endstaff;
A2   = zeros(1,n);
A2(1,[n1, 3*n1]) = [1 -1];

% Intermediate constraints
b_L3 = zeros(n1-1,1);
b_U3 = b_L3;
A3   = zeros(n1-1,n);
for i=2:n1
   A3(i-1,[i,i-1,2*n1-1+i,n1+i]) = [1,-1,1,-1];
end

% Precense constraints
b_L4 = demands;
b_U4 = demands;
A4   = zeros(n1,n);
for i=1:n1
   A4(i,[i,3*n1+i,4*n1+i]) = [1,-1,1];
end

% Overtime constraints
b_L5 = -inf*ones(n1,1);
b_U5 = zeros(n1,1);
A5   = zeros(n1,n);
for i=1:n1
   A5(i,[i, 4*n1+i]) = [-overtimemax,1];
end

% Arrival constraints
x_U(n1+1:2*n1,1) = maxtransferin*ones(n1,1);

% Leaving constraints
b_L6 = -inf*ones(n1,1);
b_U6 = zeros(n1,1);
A6   = zeros(n1,n);
for i=1:n1
   A6(i,[i, 2*n1+i]) = [-maxtransferout,1];
end

% Merging constraints

A   = [A1;A2;A3;A4;A5;A6];
b_L = [b_L1;b_L2;b_L3;b_L4;b_L5;b_L6];
b_U = [b_U1;b_U2;b_U3;b_U4;b_U5;b_U6];

c    = [zeros(n1,1); transferin*ones(n1,1); transferout*ones(n1,1); ...
        staffingdevcost*ones(2*n1,1)];

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Planning the Personnel at a Construction Site', [], [], IntVars);

% MODIFICATION LOG
%
% 051205 med   Created.