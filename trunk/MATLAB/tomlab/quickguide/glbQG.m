% glbQG is a small example problem for defining and solving
% unconstrained global programming problems using the TOMLAB format.

Name  = 'Shekel 5';
x_L   = [ 0  0  0  0]';  % Lower bounds for x.
x_U   = [10 10 10 10]';  % Upper bounds for x.
x_0   = [-3.0144 -2.4794 -3.1584 -3.1790]; % Most often not used.
x_opt = [];
f_opt = -10.1531996790582;
f_Low = -20;             % Lower bound on function.
x_min = [ 0  0  0  0]; % For plotting
x_max = [10 10 10 10]; % For plotting

Prob  = glcAssign('glbQG_f', x_L, x_U, Name, [], [], [], ... 
                   [], [], [], x_0, ...
                   [], [], [], [], ...
                   f_Low, x_min, x_max, f_opt, x_opt);

Prob.optParam.MaxFunc = 1500;              
              
Result1 = tomRun('glbFast', Prob, 1);  % Global solver            
Result2 = tomRun('conSolve', Prob, 2); % Local solver, starting from Prob.x_0
% Also possible to use a mixed-integer global solver
Result = tomRun('glcDirect', Prob, 1); 

% Result = tomRun('glbDirect', Prob, 1);
% Result = tomRun('glcDirect', Prob, 1);
% Result = tomRun('glbSolve', Prob, 1);
% Result = tomRun('glcSolve', Prob, 1);
% Result = tomRun('glcFast', Prob, 1);
% Result = tomRun('lgo', Prob, 1);
% Result = tomRun('oqnlp', Prob, 1);
