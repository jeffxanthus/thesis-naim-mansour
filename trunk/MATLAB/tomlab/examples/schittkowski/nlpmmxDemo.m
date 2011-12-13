% Demonstration of maximum norm minimization using nlpmmx
%
% Example taken from NLPMMX_demo.for by K. Schittkowski

Name = 'nlpmmxDemo (based on NLPMMX_demo.for)';

x_0 = [1.0 2.0];
x_L = [-1e5 -1e5];
x_U = [1e5 1e5];

c_L = [0.0];
c_U = [Inf];

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

Prob = clsAssign('nlpmmxDemo_r', 'nlpmmxDemo_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, 'nlpmmxDemo_c', 'nlpmmxDemo_dc', ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt, ...
                 IntVars, VarWeight, fIP, xIP);

Prob.PriLev = 2;
Prob.PriLevOpt = 4;
Prob.NLPMMX.acc = 1.0e-13;
Prob.NLPMMX.accqp = 1.0e-13;
Prob.NLPMMX.ressize = 0.0;
Prob.NLPMMX.maxfun = 20;
Prob.NLPMMX.maxit = 100;
Prob.NLPMMX.maxnm = 0;
Prob.NLPMMX.tolnm = 0.0;
Prob.NLPMMX.PrintFile = 'NLPMMX_demo.txt';
Prob.USER.P = 1;

Result = tomRun('nlpmmx',Prob);

