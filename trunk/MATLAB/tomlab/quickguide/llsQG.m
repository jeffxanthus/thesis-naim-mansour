% llsQG is a small example problem for defining and solving
% linear least squares using the TOMLAB format.
Name='LSSOL test example';           % Problem name, not required.
n = 9;                              
x_L = [-2 -2 -inf, -2*ones(1,6)]';   % Lower bounds on x
x_U = 2*ones(9,1);                   % Upper bounds on x

% Matrix defining linear constraints
A   = [ ones(1,8) 4; 1:4,-2,1 1 1 1; 1 -1 1 -1, ones(1,5)]; 
b_L = [2    -inf -4]';    % Lower bounds on the linear inequalities
b_U = [inf    -2 -2]';    % Upper bounds on the linear inequalities

% Vector m x 1 with observations in objective ||Cx -y(t)||
y = ones(10,1);  

% Matrix m x n in objective ||Cx -y(t)||
C = [ ones(1,n); 1 2 1 1 1 1 2 0 0; 1 1 3 1 1 1 -1 -1 -3; ...
      1 1 1 4 1 1 1 1 1;1 1 1 3 1 1 1 1 1;1 1 2 1 1 0 0 0 -1; ...
      1 1 1 1 0 1 1 1 1;1 1 1 0 1 1 1 1 1;1 1 0 1 1 1 2 2 3; ...
      1 0 1 1 1 1 0 2 2];

% Starting point.
x_0 = 1./[1:n]';

% x_min and x_max are only needed if doing plots.
x_min = -ones(n,1);
x_max =  ones(n,1);                          

% x_opt estimate.
x_opt = [2 1.57195927 -1.44540327 -0.03700275 0.54668583 0.17512363 ...
          -1.65670447 -0.39474418  0.31002899]; 
f_opt = 0.1390587318; % Estimated optimum.


% See 'help llsAssign' for more information.
Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
                            [], [], [], ...
                            A, b_L, b_U, ... 
                            x_min, x_max, f_opt, x_opt);
                        
Result = tomRun('clsSolve', Prob, 1);
%Result = tomRun('nlssol', Prob, 1);
%Result = tomRun('snopt', Prob, 1);
%Result = tomRun('lssol', Prob, 1);