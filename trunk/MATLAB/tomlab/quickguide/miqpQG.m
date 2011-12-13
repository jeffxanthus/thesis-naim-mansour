% miqpQG is a small example problem for defining and solving
% mixed-integer quadratic programming problems using the TOMLAB format.

c    = [-6 0]';
Name = 'XP Ref Manual MIQP';
F    = [4 -2;-2 4];
A    = [1 1];
b_L  = -Inf;
b_U  = 1.9;
x_L  = [0 0]';
x_U  = [Inf Inf]';

% Defining first variable as an integer
IntVars   = 1;

% Assign routine for defining a MIQP problem.
Prob = miqpAssign(F, c, A, b_L, b_U, x_L, x_U, [], ...
           IntVars, [], [], [], Name, [], []);

% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.

Result = tomRun('cplex', Prob, 1);
%Result = tomRun('oqnlp', Prob, 1);
%Result = tomRun('miqpBB', Prob, 1);
%Result = tomRun('xpress-mp', Prob, 1);
%Result = tomRun('minlpBB', Prob, 1);