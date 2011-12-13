% binbinQG is a small example problem for defining and solving
% combinatorial binary programming problems using the TOMLAB format.

Name  = 'binbinQG';

% The first 4 are for binary variables:
% b1 b2 b3 b4
% The rest are the combinations:
% b1*b2, b1*b3, b1*b4, b2*b3, b2*b4, b3*b4

% Coefficients in linear objective function
c     = [zeros(1,4) 1 2 3 2 4 3]';

% At least 3 binary variables are active
A     = [1 1 1 1 0 0 0 0 0 0];
b_U   = inf;
b_L   = 3;
x_L   = zeros(10,1);
x_U   = ones(10,1);

IntVars = ones(10,1);

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], Name, [], [], IntVars);

% Give the indices for the combinations and variables combined
Prob = binbin2lin(Prob, (5:10)', [1;1;1;2;2;3], [2;3;4;3;4;4]);

Result = tomRun('cplex', Prob, 1);