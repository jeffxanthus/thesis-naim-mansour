% Demonstration of least squares and constrained nonlinear data fittting using nlplsq
%
% Example taken from NLPLSQ_demo.for and NLPLSQ_demo_fit.for by K. Schittkowski

Name = 'nlplsqDemo (based on NLPLSQ_demo.for)';

x_0 = [-1.2 1.0];
x_L = [-1e5 -1e5];
x_U = [1e5 1e5];

c_L = [];
c_U = [];

JacPattern = []; ConsPattern = [];
y = zeros(length(x_0)); 
%y = [];
t = [];
weightType = []; weightY = [];
SepAlg = [];
fLowBnd = -Inf;
A = []; b_L = []; b_U = [];
x_min = []; x_max = [];
f_opt = []; x_opt = [];
IntVars = []; VarWeight = [];
fIP = []; xIP = [];

Prob = clsAssign('nlplsqDemo_r', 'nlplsqDemo_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, [], [], ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt, ...
                 IntVars, VarWeight, fIP, xIP);

Prob.PriLev = 2;
Prob.PriLevOpt = 4;
Prob.NLPLSQ.acc = 1.0e-14;
Prob.NLPLSQ.accqp = 1.0e-14;
Prob.NLPLSQ.ressize = 1.0e-14;
Prob.NLPLSQ.maxfun = 20;
Prob.NLPLSQ.maxit = 100;
Prob.NLPLSQ.maxnm = 0;
Prob.NLPLSQ.tolnm = 0.0;
Prob.NLPLSQ.PrintFile = 'NLPLSQ_demo.txt';
Prob.USER.P = 1;

Result = tomRun('nlplsq',Prob);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Name = 'nlplsqDemo (based on NLPLSQ_demo_fit.for)';

x_0 = [0.25 0.39 0.415 0.39];
x_L = [0 0 0 0];
x_U = [1e5 1e5 1e5 1e5];

c_L = [0 0];
c_U = [0 0];

JacPattern = []; ConsPattern = [];
t = [0.0625 0.0714 0.0823 0.1000 0.1250 0.1670 ...
     0.2500 0.5000 1.0000 2.0000 4.0000];
y = [0.0246 0.0235 0.0323 0.0342 0.0456 0.0627 ...
     0.0844 0.1600 0.1735 0.1947 0.1957D0];
weightType = []; weightY = [];
SepAlg = [];
fLowBnd = -Inf;
A = []; b_L = []; b_U = [];
x_min = []; x_max = [];
f_opt = []; x_opt = [];
IntVars = []; VarWeight = [];
fIP = []; xIP = [];

Prob = clsAssign('nlplsqDemo_r', 'nlplsqDemo_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, 'nlplsqDemo_c', 'nlplsqDemo_dc', ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt, ...
                 IntVars, VarWeight, fIP, xIP);

Prob.PriLev = 2;
Prob.PriLevOpt = 4;
Prob.NLPLSQ.acc = 1.0e-14;
Prob.NLPLSQ.accqp = 1.0e-14;
Prob.NLPLSQ.ressize = 1.0e-4;
Prob.NLPLSQ.maxfun = 20;
Prob.NLPLSQ.maxit = 100;
Prob.NLPLSQ.maxnm = 0;
Prob.NLPLSQ.tolnm = 0.1;
Prob.NLPLSQ.PrintFile = 'NLPLSQ_demo_fit.txt';
Prob.USER.P = 2;

Result = tomRun('nlplsq',Prob);
