% lpconQG is a small example problem for defining and solving linear
% nonlinearly constrained programming problems using the TOMLAB format.

Name = 'Linear constrained problem';
A = [1  0 -1  0
     0  1  0 -1
     0  0  1 -1];
b_L = [1 1 0]';
b_U = [];
c_L = -1;
c_U = -1;

x_0 = [1 1 1 1]';      % Starting value (not used)

x_L = [-Inf -Inf -Inf 1]';	% Lower bounds for x
x_U = 100*ones(4,1);	    % Upper bounds for x

% Objective f = 3*x(1)+2*x(2), assign as linear vector

d = [3 2 0 0]';

Prob = lpconAssign(d, x_L, x_U, Name, x_0, A, b_L, b_U,...
                   'lpconQG_c', 'lpconQG_dc', 'lpconQG_d2c', [], c_L, c_U);

% Run a global solver for 10000 function evaluations.                   
Result = tomRun('glcFast', Prob, 1);               

% Assign new starting point.
Prob.x_0 = Result.x_k(:,1);

% Run SNOPT as a local solver.

Result2 = tomRun('snopt', Prob, 1);

% Try KNITRO as a local solver, ALG = 3
% Prob.KNITRO.options.ALG = 3;
% Result3 = tomRun('knitro', Prob, 1);