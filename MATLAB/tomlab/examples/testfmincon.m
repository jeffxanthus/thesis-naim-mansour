function testfmincon

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
fprintf('Test problem from OPT TB 2.0 manual page 4.37-4.38.');
fprintf(' Two linear inequality constraints\n');
fprintf('\n');
fprintf('Initial vector is (10,10,10).');
fprintf('\n');

format long
x_0 = [10;10;10];

A=[-1 -2 2;1 2 2];
b=[0;72];

if exist('optimset')
   options=optimset;
else
   options=[];
end
options.Display='iter';
options.Diagnostics='on';
% If LargeScale = 'off' npsol will be used, if 'on' snopt (assuming Tomlab /SOL)
options.LargeScale  = 'off';  
 
[x, f_k, ExitFlag, Output, Lambda] = fmincon('fmincon_f', x_0, A, b, ...
    [], [], [], [], [], options);

fprintf('\n');
xprinte(A*x-b,            'Constraint Ax-b   ');
xprinte(Lambda.lower,     'Lambda.lower:     ');
xprinte(Lambda.upper,     'Lambda.upper:     ');
xprinte(Lambda.eqlin,     'Lambda.eqlin:     ');
xprinte(Lambda.ineqlin,   'Lambda.ineqlin:   ');
xprinte(Lambda.eqnonlin,  'Lambda.eqnonlin:  ');
xprinte(Lambda.ineqnonlin,'Lambda.ineqnonlin:');
xprinte(x,                'x:                ');
format compact
disp('Output Structure')
disp(Output)
fprintf('Function value %40.20f. ExitFlag %d\n',f_k,ExitFlag);
pause(3)

fprintf('Test problem glc1 - Floudas-Pardalos 3.3.\n');
fprintf('Four linear inequality constraints\n');
fprintf('Two quadratic (nonlinear) inequality constraints\n');
fprintf('Estimate derivatives numerically\n');
fprintf('\n');

% This example is number 16 in glc_prob.m
x_L   = [ 0  0 1 0 1 0]';  % Lower bounds on x

A =     [ 1 -3 0 0 0 0     
         -1  1 0 0 0 0
          1  1 0 0 0 0];   % Linear equations
b_L   = [-inf -inf  2 ]';  % Upper bounds for linear equations
b_U   = [  2    2   6 ]';  % Lower bounds for linear equations

%
% Original problem has no bounds on x(1) and x(2)
%
% x(1) >= 0, x(2) >= 0 & linear equation 3:  x(1) + x(2) <= 6  ==>
%
% x(1) <= 6 and x(2) <=6. This is inserted to get a box-bounded problem
x_U   = [6 6 5 6 5 10]';  % Upper bounds after x(1),x(2) values inserted
c_L   = [4 4]';           % Lower bounds on two nonlinear constraints
c_U   = [];               % Upper bounds are infinity for nonlinear constraints

x_opt = [5 1 5 0 5 10]'; 
f_opt = -310;

x_0   = (x_L+x_U)/2;      % Start with mid point

if exist('optimset')
   options=optimset;
else
   options=[];
end
options.Display     = 'iter';
options.Diagnostics = 'on';
options.GradObj     = 'off';
options.GradConstr  = 'off';
% If LargeScale = 'off' npsol will be used, if 'on' snopt (assuming Tomlab /SOL)
options.LargeScale  = 'off';  

% Rewrite for Opt Tbx format

A = [A; -A(3,:)];
b = [b_U; -b_L(3)];

n = length(x_0);

% The following row shows the criterion used by fmincon interface 
% to select Tomlab solver (unless global otxProb.Solver.Tomlab is defined)
%
% Solver = GetSolver('con',n > 200 | Prob.LargeScale,0);
%
% n = 6 and options.LargeScale = 'off' implies Prob.LargeScale == 0

Solver = GetSolver('con',0,0);

if strcmpi(Solver,'snopt')
   % fmincon interface will choose snopt (if not using opt tbx)
   Prob = ProbDef;
   Prob.Name = 'Floudas-Pardalos 3.3';
   Prob.SOL.PrintFile = 'snopt.txt';
   Prob.SOL.SummFile  = 'snopts.txt';
   Prob.SOL.optPar(1) = 11;
   Prob.SOL.optPar(11) = 1E-7;
   Prob.SOL.optPar(30) = 10000;
   Prob.SOL.optPar(35) = 1000;
   Prob.SOL.optPar(39) = 0;
   [x, f_k, ExitFlag, Output, Lambda] = fmincon('glc4_f', x_0, A, b, ...
       [],[],x_L, x_U, 'fmincon_c', options, Prob);
elseif strcmpi(Solver,'npsol')
   % fmincon interface will choose npsol (if not using opt tbx)
   Prob = ProbDef;
   Prob.Name = 'Floudas-Pardalos 3.3';
   Prob.SOL.PrintFile = 'npsol.txt';
   Prob.SOL.SummFile  = 'npsols.txt';
   Prob.SOL.optPar(1) = 11;
   Prob.SOL.optPar(13) = 3;
   Prob.SOL.optPar(30) = 1000;
   Prob.SOL.optPar(39) = 0;
   [x, f_k, ExitFlag, Output, Lambda] = fmincon('glc4_f', x_0, A, b, ...
       [],[],x_L, x_U, 'fmincon_c', options, Prob);
else
   [x, f_k, ExitFlag, Output, Lambda] = fmincon('glc4_f', x_0, A, b, ...
       [],[],x_L, x_U, 'fmincon_c', options);
end
 

fprintf('\n');
xprinte(A*x-b,            'Constraint Ax-b   ');
xprinte(Lambda.lower,     'Lambda.lower:     ');
xprinte(Lambda.upper,     'Lambda.upper:     ');
xprinte(Lambda.eqlin,     'Lambda.eqlin:     ');
xprinte(Lambda.ineqlin,   'Lambda.ineqlin:   ');
xprinte(Lambda.eqnonlin,  'Lambda.eqnonlin:  ');
xprinte(Lambda.ineqnonlin,'Lambda.ineqnonlin:');
xprinte(x,                'x:                ');
format compact
disp('Output Structure')
disp(Output)
fprintf('Function value %40.20f. ExitFlag %d\n',f_k,ExitFlag);