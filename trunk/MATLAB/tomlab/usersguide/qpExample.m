
   Name  = 'QP Example';  % File qpExample.m
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

