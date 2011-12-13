% function Prob = schedulingnurses(demand, hoursperday, breakinterval,
% breaklength, intervallength, maxdemand, flag)
%
% Creates a TOMLAB MILP problem for scheduling nurses
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 2, 2005.   Last modified Dec 2, 2005.

function Prob = schedulingnurses(demand, hoursperday, breakinterval, breaklength,...
   intervallength, maxdemand, flag)

if nargin < 7
   error('The function requires 7 inputs');
end

if isempty(demand) | isempty(hoursperday) | isempty(breakinterval) | isempty(breaklength) ...
      | isempty(intervallength) | isempty(maxdemand) | isempty(flag)
   error('One of the inputs are empty');
end

if flag == 0

   n    = length(demand); %nurses

   % FORMULATE PROBLEM
   % All variables are binary
   IntVars   = ones(n,1);
   x_L       = zeros(n,1);
   x_U       = inf*ones(n,1);

   % Requirement constraints
   b_L = demand;
   b_U = inf*ones(n,1);
   idxvec = [1:12, 1:4];
   A   = zeros(n,n);
   for i=1:n
      idx = idxvec(i:i+4);
      idx(3) = [];
      A(i,idx) = ones(1,4);
   end
   c   = ones(n,1);
end

if flag == 1

   n1    = length(demand); %nurses
   n2    = n1;             %overtime
   n     = n1+n2;

   % FORMULATE PROBLEM
   % All variables are binary
   IntVars   = ones(n,1);
   x_L       = zeros(n,1);
   x_U       = inf*ones(n,1);

   % Max demand constraint
   b_L1 = -inf;
   b_U1 = maxdemand;
   A1   = [ones(1,n1), zeros(1,n2)];

   % Requirement constraints
   b_L2 = demand;
   b_U2 = inf*ones(n1,1);
   idxvec1 = [1:12, 1:4];
   idxvec2 = [1:12, 1:5]+n1;
   A2   = zeros(n1,n);
   for i=1:n1
      idx = idxvec1(i:i+4);
      idx(3) = [];
      A2(i,idx) = ones(1,4);
      idx = idxvec2(i+5);
      A2(i,idx) = ones(1,1);
   end
   
   % Overtime less than regular (extends these)
   A3 = zeros(n1,n);
   b_L3 = -inf*ones(n1,1);
   b_U3 = zeros(n1,1);
   for i=1:n1
      A3(i,[i,i+n1]) = [-1 1];
   end
   
   % Merge constraints
   A   = [A1;A2;A3];
   b_L = [b_L1;b_L2;b_L3];
   b_U = [b_U1;b_U2;b_U3];
   
   c   = [zeros(n1,1);ones(n2,1)];
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Scheduling Nurses', [], [], IntVars);

% MODIFICATION LOG
%
% 051202 med   Created.