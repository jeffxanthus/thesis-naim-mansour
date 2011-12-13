% Demonstration of multicriteria optimization using nlpjob
%
% Example taken from Demo.for by K. Schittkowski


% Demo.for

Name = 'nlpjobDemo (Demo.for)';

x_0 = [1 1];
x_L = [-10 -10];
x_U = [10 10];

c_L = [0 0];
c_U = [Inf Inf];

JacPattern = []; ConsPattern = [];
y = zeros(length(x_0)); t = [];
weightType = []; weightY = [];
SepAlg = [];
fLowBnd = -Inf;
A = []; b_L = []; b_U = [];
x_min = []; x_max = [];
f_opt = []; x_opt = [];
IntVars = []; VarWeight = [];
fIP = []; xIP = [];

Prob = clsAssign('nlpjobDemo_f', 'nlpjobDemo_g', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, 'nlpjobDemo_c', 'nlpjobDemo_dc', ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt, ...
                 IntVars, VarWeight, fIP, xIP);

Prob.NLPJOB.acc = 1e-10;
Prob.NLPJOB.accqp = 1e-10;
Prob.NLPJOB.maxit = 100;
Prob.NLPJOB.maxfun = 10;
Prob.NLPJOB.w = [10 10];
Prob.NLPJOB.fk = [1 -3];
Prob.NLPJOB.imin =  2;
Prob.NLPJOB.PrintFile = 'nlpjobDemo.txt';
Prob.PriLevOpt = 0;
Prob.PriLev = 0;
Prob.USER.P = 1;

for i = 0:15
   Prob.NLPJOB.model = i;
   if i == 3
      Prob.NLPJOB.w(1) = 2;
   end
   disp(['Testing nlpjob with model = ' num2str(i)]);
   Prob.NLPJOB.PrintFile = ['nlpjobDemo_Model_' num2str(i) '.txt'];
   Prob.PriLev = 0;
   Prob.PriLevOpt = 3;
   Result = tomRun('nlpjob',Prob);
end
