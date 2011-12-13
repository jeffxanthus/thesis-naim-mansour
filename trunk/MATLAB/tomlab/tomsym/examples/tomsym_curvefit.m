%% Least squares example

% Data from Yalmip example
x = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)';
a = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4+8*rand(length(a),1));
e(100:115) = 30;
y = a*x+e;

% Define the decision variable
toms 6x1 x_hat

% x_hat and the regressors a define the residuals with y
residuals = y-a*x_hat;

%The L2 problem is solved as a QP problem without any constraints.
options = struct;
options.norm = 'L2';
solution = ezsolve(residuals,[],[],options);

plot(t,subs(a*x_hat,solution),'-',t,y,'.')