% simQG is a small example problem for defining and solving simulation
% problems where the objective function and constraints are evaluated
% during the same function call.

Name = 'HS 47';
b_L  = [];
b_U  = [];
A    = [];
c_L  = [0; 0; 0];
c_U  = [0; 0; 0];
x_0  = [2; sqrt(2); -1; 2-sqrt(2); .5];
x_L  = [];
x_U  = [];

Prob = simAssign('simQG_fc', 'simQG_gdc', [], [], x_L, x_U, ...
       Name, x_0, [], A, b_L, b_U, [], c_L, c_U);
   
Result = tomRun('snopt', Prob, 1);
% Prob.KNITRO.options.HESSOPT = 6;
% Prob.KNITRO.options.ALG = 3;
% Result2 = tomRun('knitro', Prob, 1);