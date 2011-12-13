% qpconQG is a small example problem for defining and solving quadratic
% nonlinearly constrained programming problems using the TOMLAB format.

Name  = 'QP constrained problem';
x_0   = ones(10,1);
x_L   = [];
x_U   = [];
A     = ones(8,10);
for i = 1:8
    A(i,i) = 1/2;
end
b_L   = ones(8,1);
b_U   = ones(8,1);
c_L   = 4;
c_U   = 4;

% Objective f = -x'*x + sum(x), assign as quadratic/linear matrix/vector

F     = -2*speye(10);
d     = ones(10,1);

Prob = qpconAssign(F, d, x_L, x_U, Name, x_0, A, b_L, b_U,...
    'qpconQG_c', 'qpconQG_dc', 'qpconQG_d2c', [], c_L, c_U);

% Run SNOPT as a local solver
Result = tomRun('snopt', Prob, 1);
% Result2 = tomRun('minos', Prob, 1);
% Result3 = tomRun('npsol', Prob, 1);
% Result4 = tomRun('knitro', Prob, 1);