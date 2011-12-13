% gpQG is a small example problem for defining and solving
% geometric programming problems using the TOMLAB format.

nterm = [6;3];
coef = [.5e1;.5e5;.2e2;.72e5;.1e2;.144e6;.4e1;.32e2;.12e3];
A = sparse([ 1  -1  0   0  0  0 -1  0  0;...      
             0   0  1  -1  0  0  0 -1  0;...
             0   0  0   0  1 -1  0  0 -1])';
                  
Name  = 'GP Example';  % File gpQG.m

% Assign routine for defining a GP problem.
Prob = gpAssign(nterm, coef, A, Name);

% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.

Result = tomRun('GP', Prob, 1);