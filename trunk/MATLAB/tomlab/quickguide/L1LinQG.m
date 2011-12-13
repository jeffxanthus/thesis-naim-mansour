% L1LinQG is a small example problem for defining and solving
% a linearly constrained linear L1 problem using the TOMLAB format.

Name='L1LinSolve test example';      % Problem name, not required.
n = 6;                              
x_L = -10*ones(n,1);                 % Lower bounds on x
x_U =  10*ones(n,1);                 % Upper bounds on x
x_0 = (x_L + x_U) / 2;               % Starting point

C = spdiags([1 2 3 4 5 6]', 0, n, n); % C matrix
y = 1.5*ones(n,1);                    % Data vector

% Matrix defining linear constraints
A = [1 1 0 0 0 0];
b_L = 1;                  % Lower bounds on the linear inequalities
b_U = 1;                  % Upper bounds on the linear inequalities

% Defining damping matrix
Prob.LS.damp = 1;
Prob.LS.L = spdiags(ones(6,1)*0.01, 0, 6, 6);

% See 'help llsAssign' for more information.
Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
                            [], [], [], ...
                            A, b_L, b_U);
                        
Prob.SolverL1 = 'lpSimplex';
Result = tomRun('L1LinSolve', Prob, 1);
% Prob.SolverL1 = 'MINOS';
% Result = tomRun('L1LinSolve', Prob, 1);
% Prob.SolverL1 = 'CPLEX';
% Result = tomRun('L1LinSolve', Prob, 1);