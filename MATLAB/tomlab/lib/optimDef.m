% Define parameter vector optPar used as input to solvers from Optimization TB
% See Matlab file foptions.m for a description of the elements
%
% function optPar = optimDef(optParam, LineParam, Prob);
%
% INPUT:
%  optParam   Structure
%  LineParam  Structure
%
% OUTPUT:
%  optPar     OPTIONS vector with length 18
%
%  optPar(13), the number of equality constraints, is be computed elsewhere

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Sept 12, 1998.   Last modified Nov 6, 2000.

function optPar = optimDef(optParam, LineParam, Prob)

optPar=zeros(1,18);

% (1) Print level in optimization solver
optPar(1)=Prob.PriLevOpt;

% (2) Convergence tolerance in optimal solution x
% Distance between successive x. //x_k+1 - x_k//
optPar(2)=optParam.eps_x;

% (3) Convergence tolerance on F.(Default=1E-10). 
% In OPTIM TB, used for directed.derivative:  g_k^T * p_k <= eps_f 
optPar(3)=optParam.eps_dirg;

% (4) Constraint violation convergence tolerance
optPar(4)=optParam.cTol;

% (5) Optimization Algorithm. Dependent on type of problem. Set as default 0.
% NOT USED, instead Prob.Solver.Alg optPar(5)=optParam.alg;

% (6) Optimization Method. Subalgorithm selection. Problem dependent. Set to 0.
optPar(6)=Prob.Solver.Method;

% (7) Line search algorithm. 0 = quadratic, 1 = cubic
optPar(7)=LineParam.LineAlg;

% OPTIONS(8)  - Function value. (Lambda in goal attainment. )

% (9) Set to 1 if you want to check user-supplied gradients
optPar(9)=optParam.GradCheck > 0;

% OPTIONS(10) - Number of Function and Constraint Evaluations.
% OPTIONS(11) - Number of Function Gradient Evaluations.
% OPTIONS(12) - Number of Constraint Evaluations
% OPTIONS(13) - Number of equality constraints. 
% This number is implicit, i.e.
% me = OPTIONS(13)= sum(Prob.b_L==Prob.b_U) +  sum(Prob.c_L==Prob.c_U)

% me is computed elsewhere

% (14) Maximal number of iterations
optPar(14)=optParam.MaxIter;

% OPTIONS(15) - Used in goal attainment for special objectives. 

% (16) Minimum change in variables for finite difference gradients.
optPar(16)=0.01*optParam.DiffInt;

% (17) Maximum change in variables for finite difference gradients.
optPar(17)=0.1;

% (18) Initial step length. (Default 1 or less). 
optPar(18)=LineParam.InitStepLength;