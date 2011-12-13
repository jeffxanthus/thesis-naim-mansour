% nlpQG is a small example problem for defining and solving
% nonlinear programming problems using the TOMLAB format.

Name = 'RBB Problem';
x_0 = [-1.2 1]';     % Starting values for the optimization.
x_L = [-10;-10];     % Lower bounds for x.
x_U = [2;2];         % Upper bounds for x.
fLowBnd = 0;         % Lower bound on function.

c_L = -1000;         % Lower bound on nonlinear constraints.
c_U = 0;             % Upper bound on nonlinear constraints.

Prob = conAssign('rbbQG_f', 'rbbQG_g', 'rbbQG_H', [], x_L, x_U, Name, x_0,...
    [], fLowBnd, [], [], [], 'rbbQG_c', 'rbbQG_dc', 'rbbQG_d2c', [], c_L, c_U);

Prob.Warning = 0;    % Turning off warnings.            
            
Result = tomRun('ucSolve', Prob, 1);  % Ignores constraints.

% Result = tomRun('conopt', Prob, 1);
% Result = tomRun('snopt', Prob, 1);