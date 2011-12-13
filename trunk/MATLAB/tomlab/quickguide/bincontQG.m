% bincontQG is a small example problem for defining and solving
% where a subset of the continuous variables are active.

Name  = 'bincontQG';

% The first n are for binary variables:
% b1 ... bn
%
% Variables n+1 to 2*n are continuous:
% c1 ... cn
%
% Variables 2*n+1 to 3*n are combinations:
% b1*c1 ... bn*cn

% Coefficients in linear objective function
n   = 20;
cov1 = magic(20)/100;
cov1 = cov1(:,1);
F   = [zeros(2*n,3*n); zeros(n,2*n), diag(cov1)];
c   = zeros(3*n,1);

% 10 variables to be selected
% Combined combinations greater than 1000
A     = [ones(1,n), zeros(1,2*n);
    zeros(1,2*n), ones(1,n)];
b_U   = [10;inf];
b_L   = [10;1000];
x_L   = zeros(3*n,1);
x_U   = [ones(n,1); 1e4*ones(2*n,1)];

IntVars = [ones(n,1); zeros(2*n,1)];

Prob = miqpAssign(F, c, A, b_L, b_U, x_L, x_U, [], IntVars);

% Give the indices for the combinations and variables combined
Prob = bincont2lin(Prob, (2*n+1:3*n)', (1:n)', (n+1:2*n)');

Result = tomRun('cplex', Prob, 1);