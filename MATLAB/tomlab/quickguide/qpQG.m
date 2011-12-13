% qpQG is a small example problem for defining and solving
% quadratic programming problems using the TOMLAB format.

Name  = 'QP Example';
F     = [ 8   1        % Matrix F in 1/2 * x' * F * x + c' * x
          1   8 ];
c     = [ 3  -4 ]';    % Vector c in 1/2 * x' * F * x + c' * x
A     = [ 1   1        % Constraint matrix
          1  -1 ];
b_L   = [-inf  0  ]';  % Lower bounds on the linear constraints
b_U   = [  5   0  ]';  % Upper bounds on the linear constraints
x_L   = [  0   0  ]';  % Lower bounds on the variables
x_U   = [ inf inf ]';  % Upper bounds on the variables
x_0   = [  0   1  ]';  % Starting point
x_min = [-1 -1 ];      % Plot region lower bound parameters
x_max = [ 6  6 ];      % Plot region upper bound parameters

% Assign routine for defining a QP problem.
Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name,...
       [], [], [], x_min, x_max);

% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.

Result = tomRun('qpSolve', Prob, 1);
%Result = tomRun('snopt', Prob, 1);
%Result = tomRun('sqopt', Prob, 1);
%Result = tomRun('cplex', Prob, 1);
%Result = tomRun('knitro', Prob, 1);
%Result = tomRun('conopt', Prob, 1);