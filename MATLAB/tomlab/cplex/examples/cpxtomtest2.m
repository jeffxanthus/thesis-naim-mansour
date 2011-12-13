% Test of TOMLAB / CPLEX call using TOMLAB input format
%
% Simple small problem solved as IP using TOMLAB (TQ) format

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2006 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written July 28, 2001.       Last modified Mar 22, 2005.

format compact
fprintf('=====================================================\n');
fprintf('Run very simple IP defined in the TOMLAB format\n');
fprintf('=====================================================\n');

A=[8 5 
  -1 2];

disp('The matrix A')
A

disp('The right hand side b')
b=[31 6]';       % right-hand side

b

disp('The objective function coefficients c')
c=[-30  -18 ]'; % cost vector 

c

disp('Set lower bounds for x as zero and upper bounds infinite')

[m,n] = size(A);

x_L = zeros(n,1);       % NOTE! Must be set as 0, otherwise assumed -Inf
x_U = [];               % Default set as infinite


% Use the TOMLAB (TQ) format
%
% Call the TOMLAB mipAssign routine that creates a structure with all
% problem information.

Prob = mipAssign(c, A, [], b, x_L, x_U); 

% Call the CPLEX solver

fprintf('-----------------------------------------------------\n');
fprintf('\n\nIf not setting any variables to be integers,\n');
fprintf('the problem is solved as a linear programming (LP) problem\n\n');
fprintf('-----------------------------------------------------\n');

Result = tomRun('cplex',Prob,2);

pause
fprintf('Set all variables as integer valued, and rerun\n');

% All variables should be integer
IntVars=[1:n];

fprintf('Define the name of the problem as "Simple MIP example":\n');
Name = 'Simple MIP example';

Prob = mipAssign(c, A, [], b, x_L, x_U, [], Name, [], [], IntVars);


Result = tomRun('cplex',Prob,2);
pause

fprintf('===============================================================\n');
fprintf('Again run the very simple IP defined in the TOMLAB format\n');
fprintf('\nThis time add slack variables!\n');
fprintf('\nAnd define Name and integer values directly in mipAssign\n');
fprintf('===============================================================\n');

A=[8 5 1 0
  -1 2 0 1 ];

disp('The matrix A')
A

[m,n] = size(A);

disp('The objective function coefficients c')
c=[-30  -18  0 0]'; % cost vector 
c

% All variables should be integer
IntVars=[1:n];    

x_L = zeros(n,1);      

Name = 'Simple MIP example';

% Initialize all additional variables as empty (default) and explain them

b_U   = b;      % Upper bounds set as b
b_L   = b;      % Lower bounds equal to upper bounds, now when slacks
                % have been added
x_min = [];     % Bounds for plotting, normally not used for MIP
x_max = [];     % Bounds for plotting, normally not used for MIP
x_0   = [];     % Initial value not needed

f_opt = [];     % The function value at the optimal point
                % Used in the printout. Assumed not known
x_opt = [];     % The optimal integer solution is assumed not to be known

f_Low = [];     % f_Low <= f_optimal must hold
                % Set as -realmax, lowest possible number

fIP   = [];     % Do not use any prior knowledge about function values
                % at some point where the solution is integer
xIP   = [];     % Do not use any prior knowledge about integer points

setupFile = []; % Just define the Prob structure, not any permanent Init File,
                % i.e. do not create a file according to the
                % TOMLAB Init File (IF) format

nProblem  = []; % Number of problems in the IF format is 0

VarWeight = []; % No variable priorities, largest fractional part will be used
Knapsack  = []; % No knapsack problem.

% Use the TOMLAB Quick (TQ) format with all input arguments

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, Knapsack, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);

% Call the CPLEX solver using the TOMLAB driver, add printout

Result = tomRun('cplex',Prob, 2);

% MODIFICATION LOG:
%
% 050113 med Modification log added
% 050209 med Changed name from tomtest2 to cpxtomtest2
