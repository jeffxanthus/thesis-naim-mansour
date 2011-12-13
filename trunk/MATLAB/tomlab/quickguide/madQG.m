% madQG are two examples for defining and solving nonlinear
% programming problems using TOMLAB /MAD

Name = 'RBB Problem';
x_0 = [-1.2 1]';     % Starting values for the optimization
x_L = [-10;-10];     % Lower bounds for x.
x_U = [2;2];         % Upper bounds for x.
fLowBnd = 0;         % Lower bound on function.

c_L = -1000;         % Lower bound on nonlinear constraints.
c_U = 0;             % Upper bound on nonlinear constraints.

Prob1 = conAssign('rbbQG_f', [], [], [], x_L, x_U, Name, x_0,...
                [], fLowBnd, [], [], [], 'rbbQG_c', [], [], [], c_L, c_U);

Prob2 = conAssign('rbbQG_f', 'rbbQG_g', [], [], x_L, x_U, Name, x_0,...
                [], fLowBnd, [], [], [], 'rbbQG_c', 'rbbQG_dc', [], [], c_L, c_U);
                             
Prob1.Warning = 0;    % Turning off warnings.            
Prob2.Warning = 0;    % Turning off warnings.            
            
madinitglobals;
Prob1.ADObj  = 1; % Gradient calculated
Prob1.ADCons = 1; % Jacobian calculated
Result1 = tomRun('snopt', Prob1, 1);  % Only uses first order information.

madinitglobals;
Prob2.CONOPT.options.LS2PTJ = 0;
Prob2.ADObj  = -1; % Hessian calculated
Prob2.ADCons = -1; % Lagrangian function for the nonlinear constraints.
Result2 = tomRun('conopt', Prob2, 1);  % Uses second order information.