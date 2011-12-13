mipExample;

nProblem  = 7;  % Use the same problem number as in mip_prob.m
fIP       = []; % Do not use any prior knowledge
xIP       = []; % Do not use any prior knowledge
setupFile = []; % Just define the Prob structure, not any permanent setup file
x_opt     = []; % The optimal integer solution is not known
VarWeight = []; % No variable priorities, largest fractional part will be used
KNAPSACK  = 0;  % First run without the knapsack heuristic

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);

Prob.Solver.Alg       = 2;    % Depth First, then Breadth (Default Depth First)
Prob.optParam.MaxIter = 5000; % Must increase iterations from default 500
Prob.optParam.IterPrint = 0;
Prob.PriLev = 1;
Result                = tomRun('mipSolve', Prob, 0);

% ------------------------------------------------
% Add priorities on the variables
% ------------------------------------------------
VarWeight = c;   
% NOTE. Prob is already defined,  could skip mipAssign, directly set:
% Prob.MIP.VarWeight=c;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);


Prob.Solver.Alg       = 2;             % Depth First, then Breadth search
Prob.optParam.MaxIter = 5000;          % Must increase number of iterations
Prob.optParam.IterPrint = 0;
Prob.PriLev = 1;
Result                = tomRun('mipSolve', Prob, 0);

% ----------------------------------------------
% Use the knapsack heuristic, but not priorities
% ----------------------------------------------
KNAPSACK  = 1; VarWeight = []; 
% NOTE. Prob is already defined,  could skip mipAssign, directly set:
% Prob.MIP.KNAPSACK=1;
% Prob.MIP.VarWeight=[];

Prob      = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, ...
                      nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                      f_Low, x_min, x_max, f_opt, x_opt);

Prob.Solver.Alg = 2;              % Depth First, then Breadth search
Prob.optParam.IterPrint = 0;
Prob.PriLev = 1;
Result                = tomRun('mipSolve', Prob, 0);

% --------------------------------------------------
% Now use both the knapsack heuristic and priorities
% --------------------------------------------------
VarWeight = c; KNAPSACK  = 1;
% NOTE. Prob is already defined,  could skip mipAssign, directly set:
% Prob.MIP.KNAPSACK=1;
% Prob.MIP.VarWeight=c;

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, nProblem,...
                 IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                 f_Low, x_min, x_max, f_opt, x_opt);

Prob.Solver.Alg = 2;              % Depth First, then Breadth search
Prob.optParam.IterPrint = 0;
Prob.PriLev = 1;
Result                = tomRun('mipSolve', Prob, 0);
