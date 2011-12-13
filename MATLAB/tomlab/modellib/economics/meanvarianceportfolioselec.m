% function Prob = meanvarianceportfolioselec(estreturn, covmat, target)
%
% Creates a TOMLAB MIQP problem for mean variance portfolio selection
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 2, 2005.   Last modified Dec 2, 2005.

function Prob = meanvarianceportfolioselec(estreturn, covmat, target, maxassets)

if nargin < 3
   error('The function requires 3 inputs');
end

if isempty(estreturn) | isempty(covmat) | isempty(target)
   error('One of the inputs are empty');
end

if nargin < 4
   n   = length(estreturn);  %returns

   % FORMULATE PROBLEM
   % No variables are integer, they are fractions
   IntVars   = zeros(n,1);
   x_L       = zeros(n,1);
   x_U       = ones(n,1);

   % Cost constraints
   b_L = [1;target];
   b_U = [1;inf];
   A   = [ones(1,n);estreturn'];

   c   = zeros(n,1);
   F   = covmat;
else
   if isempty(maxassets)
      error('One of the inputs are empty');
   end
   n1 = length(estreturn);
   n  = n1*2;
   % FORMULATE PROBLEM
   % No variables are integer, they are fractions
   IntVars   = [zeros(n1,1);ones(n1,1)];
   x_L       = zeros(n,1);
   x_U       = ones(n,1);

   % Cost constraints
   b_L = [1;target;-inf;-inf*ones(n1,1)];
   b_U = [1;inf;maxassets;zeros(n1,1)];
   A   = [ones(1,n1), zeros(1,n1) ;estreturn', zeros(1,n1);...
          zeros(1,n1), ones(1,n1); eye(n1), -eye(n1)];
   c   = zeros(n,1);
   F   = [covmat, zeros(n1,n1); zeros(n1,2*n1)];
end

Prob = miqpAssign(F, c, A, b_L, b_U, x_L, x_U, [], ...
   IntVars, [], [], [], 'Mean Variance Portfolio Selection');

% MODIFICATION LOG
%
% 051202 med   Created.