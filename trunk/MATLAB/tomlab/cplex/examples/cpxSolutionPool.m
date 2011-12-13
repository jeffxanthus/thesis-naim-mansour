% Test of TOMLAB /CPLEX solution pool capabilities.
%
% function [x,f] = cpxSolutionPool
%
% Exemplifies the use of the solution pool features in TOMLAB /CPLEX and
% the parameters associated with this, SOLNPOOLCAPACITY, SOLNPOOLGAP,
% SOLNPOOLINTENSITY and SOLNPOOLREPLACE.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 11.2.0$
% Written Nov 26, 2008.   Last modified Nov 26, 2008.

function [x,f] = cpxSolutionPool

% Create a test problem from the TOMLAB sets
Prob = probInit('mip_prob', 15);

% Set the solution capacity to 5. This means that 5+1 solutions will be
% available, where the first one is the best.
Prob.MIP.cpxControl.SOLNPOOLCAPACITY = 5;

% Solve the problem with CPLEX using a print level of 0
Result = tomRun('cplex', Prob, 0);
n = length(Result.f_k);

disp(sprintf('With a capacity of 5, %d solutions were generated.',n));
disp(' ');

% Now set the gap to 0.05 to avoid very bad solutions, i.e. that are more
% than 5 % from the best solution
Prob.MIP.cpxControl.SOLNPOOLGAP = 0.05;
Result = tomRun('cplex', Prob, 0);

maxgap = (min(Result.f_k) - max(Result.f_k))/min(Result.f_k);
disp(sprintf('With a solution gap of 0.05, the worst deviation is %d.',maxgap));
disp(' ');

% Set the pool intensity to very aggressive
Prob.MIP.cpxControl.SOLNPOOLINTENSITY = 4;
Prob.MIP.cpxControl.SOLNPOOLCAPACITY = 50;
Result = tomRun('cplex', Prob, 0);
n = length(Result.f_k);
disp(sprintf('With aggressive enumeration, %d solutions were found.',n));
disp(' ');

% Try to build a diverse set
Prob.MIP.cpxControl.SOLNPOOLREPLACE = 3;
Prob.MIP.cpxControl.SOLNPOOLGAP = 1e75;
Result = tomRun('cplex', Prob, 0);

maxgap = (min(Result.f_k) - max(Result.f_k))/min(Result.f_k);
disp(sprintf('Building a diverse set generates a gap of %d.',maxgap));

% The x and f are generated column-wise, i.e. each column is one solution.
% The first column is always the best solution.
x = Result.x_k(:,1);
f = Result.f_k(1);

% MODIFICATION LOG:
%
% 081126 med Written