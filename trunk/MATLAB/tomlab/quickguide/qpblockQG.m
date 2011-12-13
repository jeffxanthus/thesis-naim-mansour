% qpblockQG is a small example problem for defining and solving
% quadratic programming on a block structure using the TOMLAB format.
% See help qpblockAssign for more information.
Name  = 'QP Block Example';

switch 1
    case 1
        F     = [ 8   0
                  0   8 ];
        Fb(1).out = ones(3,2);
        Fb(1).inn = ones(3,3);
        Fb(2).out = ones(4,2);
        Fb(2).inn = ones(4,4);
        d     = [ 3  -4 ]';
    case 2
        F     = [ 8   0
                  0   8 ];
        Fb(1).out = ones(3,2);
        Fb(2).out = ones(4,2);
        d     = [ 3  -4 ]';
    case 3
        F     = [ 8   8]';
        Fb(1).out = ones(3,2);
        Fb(1).inn = ones(3,3);
        Fb(2).out = ones(4,2);
        Fb(2).inn = ones(4,4);
        d     = [ 3  -4 ]';
    case 4
        F     = [ 8   8]';
        Fb(1).out = ones(3,2);
        Fb(2).out = ones(4,2);
        d     = [ 3  -4 ]';
end

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
Prob = qpblockAssign(F, Fb, d, x_L, x_U, Name, x_0, A, b_L, b_U);

% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.

Result = tomRun('snopt', Prob, 1);
% Result = tomRun('knitro', Prob, 1);
% Result = tomRun('conopt', Prob, 1);