   Name  = 'lptest';
   c     = [-7 -5]';   % Coefficients in linear objective function 
   A     = [ 1  2
             4  1 ];   % Matrix defining linear constraints
   b_U   = [ 6 12 ]';  % Upper bounds on the linear inequalities
   x_L   = [ 0  0 ]';  % Lower bounds on x

   % x_min and x_max are only needed if doing plots
   x_min = [ 0  0 ]';
   x_max = [10 10 ]';

   % b_L, x_U and x_0 have default values and need not be defined.
   % It is possible to call lpAssign with empty [] arguments instead
   b_L   = [-inf -inf]';
   x_U   = [];
   x_0   = [];
