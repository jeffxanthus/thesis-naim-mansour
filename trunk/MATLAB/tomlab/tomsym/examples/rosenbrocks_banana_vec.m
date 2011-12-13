% rosenbrocks_banana - TomSym NLP demonstration
% In this example, x is a 2x1 vector. This is NOT a good way of using
% tomSym, because derivatives of functions involving x(1) and x(2) will be
% much more complicated than if x1 and x2 had been separate variables.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc.

toms 2x1 x
a = 100;

% Objective function
f = a*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Constraints
c   = -x(1)^2 - x(2);
con = { -1000 <= c <= 0; -10 <= x <= 2};

% Initial conditions
x0 = struct('x',[-1.2; 1]);
options = struct;
options.name = 'Rosenbrocks banana';
solution = ezsolve(f,con,x0,options);

disp('x =');
disp(solution.x);