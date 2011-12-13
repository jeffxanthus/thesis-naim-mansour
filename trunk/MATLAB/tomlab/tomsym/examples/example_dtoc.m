% Example of discrete-time optimal control using tomSym.
%
% This is an adaptation of the YALMIP example (Standatd MPC) modified to
% use tomSym.

% Model data
A = [2 -1;1 0];
B = [1;0];
nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q = eye(2);
R = 1.0001; % .0001 is tie-breaker for u4, to avoid randomness in solution.
N = 4;

% Initial state
x0 = [2;1];

constraints = {};
u = cell(1,N);
objective = 0;
x = x0;
for k = 1:N
    u{k} = tom(['u' num2str(k)], nu, 1);
    x = A*x + B*u{k};
    %objective = objective + sum((Q*x).^2) + sum((R*u{k}).^2); %QP
    objective = objective + norm(Q*x,1) + norm(R*u{k},1);      %LP
    constraints = {constraints{:}, -1 <= u{k}<= 1, -5<=x<=5};
end

solution = ezsolve(objective,constraints);