% TOMLAB DFNLP Solver
%
% function Result = dfnlpTL(Prob)
%
% INPUT:
%
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% c_L, c_U  Bounds on nonlinear constraints.
% A         Linear constraint matrix.
% PriLevOpt Print level in solver.
% WarmStart If true, use warm start, otherwise cold start.
%
% -----------------------------------------------
% Fields used in Prob.DFNLP:
% -----------------------------------------------
% maxfun    The integer variable defines an upper bound for the number
%           of function calls during the line search.
%
% maxit     Maximum number of outer iterations, where one iteration
%           corresponds to one formulation and solution of the quadratic
%           programming subproblem, or, alternatively, one evaluation of
%           gradients.
%
% acc       The user has to specify the desired final accuracy (e.g. 1.0e-7).
%           The termination accuracy should not be smaller
%           than the accuracy by which gradients are computed.
%
% ressiz    The user must indicate a guess for the approximate size of the least squares
%           residual, i.e. a low positive real number if the residual is supposed to be small,
%           and a large one in the order of 1 if the residual is supposed to be large. If
%           model is not equal to 2, ressiz must not be set by the user.
%
% model     1-4. See the user's guide for more information.
%           1: L1 DATA FITTING.
%           2: L2 - OR LEAST SQUARES DATA FITTING.
%           3: MAXIMUM-NORM DATA FITTING.
%           4: MAXIMUM FUNCTION.
%
% PrintFile Name of DFNLP Print file. Amount and type of printing determined
%           by PriLevOpt.
%
% -----------------------------------------------------------------------
%
% OUTPUT:
%
% Result   Structure with results (see ResultDef.m):
%
% x_k      Solution vector.
% x_0      Initial solution vector.
%
% f_k      Function value at optimum.
% g_k      Gradient of the objective function.
%
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
%
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
%
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
%
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
%
% ExitFlag Exit status.
% ExitText Exit text from solver.
%
% Inform   DFNLP information parameter.
%
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% Solver   Name of the solver (DFNLP).
% SolverAlgorithm  Description of the solver.
%
% DFNLP.act  The logical array indicates constraints, which DFNLP considers to be
%            active at the last computed iterate.
%
% DFNLP.u    Contains the multipliers with respect to the actual iterate
%            stored in x_k. See user's guide.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release : 4.7.0$
% Written Jun 1, 2003.    Last modified Dec 2, 2004.

function Result = dfnlpTL(Prob)

%#function dffunc dfgrad

if nargin < 1, error('dfnlpTL needs the Prob structure as input');end

global ksNLP_x ksNLP_xJ % Max number of variables/constraints/resids to print

ksNLP_xJ = [];

if isempty(Prob.LS)
   fprintf('Prob.LS not defined\n')
   error('No least squares or approximation problem')
end
%if any(Prob.probType ==[2 7 8 11])
%   error('dfnlp cannot solve LP, MILP, QP and MIQP problems');
%end

Prob.solvType = 6;          % constrained LS solver, cls

Prob = iniSolve(Prob,3,1,1);% Init globals

Prob = ksDef(Prob);         % Klaus Schittkowski problem definitions

if ~isfield(Prob,'DFNLP')  % Make sure DFNLP field exists
   Prob.DFNLP = [];
end

Result=ResultDef(Prob);     % Define the Result struct
Result.Solver='DFNLP';     % Solver name
Result.SolverAlgorithm='Dense SQP';

PriLev=Prob.PriLevOpt;      % Printing level in solver

%
% Define lower and upper bound arrays for NPSOL
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints
%

% Initial checks on the inputs

BIG=1E12;

% Use this function, even though we later separate variable bounds 
% from the rest of the constraints
[bl, bu, n, m1, m2] = defblbu(Prob, BIG,1);

% Initial checks on the inputs

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

m = m1+m2;

xl = bl(1:n); % Prob.x_L
xu = bu(1:n); % Prob.x_U

me     = Prob.mEQ;

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( xl,min(xu,x_0(:) ) ); 

% [f, g ] = nlresid(x_0, Prob);
r_0 = nlresid(x_0, Prob);
L   = length(r_0);  % Number of objective functions (or residuals)

% Define default solver options.
maxfun    = DefPar(Prob.DFNLP,'maxfun',400);
maxit     = DefPar(Prob.DFNLP,'maxit',2000);
acc       = DefPar(Prob.DFNLP,'acc',1e-6);
ressiz    = DefPar(Prob.DFNLP,'ressiz',1e-14);
PrintFile = DefPar(Prob.DFNLP,'PrintFile','dfnlp.txt');
options   = [maxfun, maxit, acc, ressiz, eps];

% Default objective functions
model     = DefPar(Prob.DFNLP,'model',1);

[x,Inform,r_k,f_k,dres,g_k,dg_k,u,act] = dfnlp(m,n,L,me,x_0,xl,xu,PriLev,PrintFile,options,model,Prob);

Result.f_k  = 0.5*f_k;
Result.r_k  = r_k; 

Result.g_k  = dres;

ksNLP_x = []; %Reset global

c_k   = zeros(m2,1);
dc_k  = zeros(m2,n);
cixEQ = Prob.cixEQ;
if ~isempty(cixEQ)
   c_k(cixEQ)    = g_k(cixEQ) + Prob.c_L(cixEQ);
   dc_k(cixEQ,:) = dg_k(cixEQ,:);
end
cixLow = Prob.cixLow;
if ~isempty(cixLow)
   c_k(cixLow)    = g_k(cixLow) + Prob.c_L(cixLow);
   dc_k(cixLow,:) = dg_k(cixLow,:);
end
cixUpp = Prob.cixUpp;
if ~isempty(cixUpp)
   c_k(cixUpp)    = Prob.c_U(cixUpp) - g_k(cixUpp);
   dc_k(cixUpp,:) = -dg_k(cixUpp,:);
end

if m1>0, Result.Ax  = Prob.A*x; else Result.Ax = [];  end

Result.c_k  = c_k;
Result.cJac = dc_k;
Result.x_k  = x;

optParam = Prob.optParam;

Result = StateDef(Result, x(1:n), Result.Ax, Result.c_k, ...
   optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);

Result.DFNLP.u = u;
Result.DFNLP.act = act;

ExitFlag = 99;
switch(Inform)
   case 0
      ExitFlag = 0;
      ExitText = 'The optimality conditions are satisfied.';
   case 1 % iterlimit
      ExitFlag = 1;
      ExitText = 'The algorithm has been stopped after MAXIT iterations.';
   case 2
      ExitText = 'The algorithm computed an uphill search direction.';
   case 3
      ExitText = 'Underflow occurred.';
   case 4
      ExitText = 'More than maxfun function evaluations are required during line search.';
   case 5
      ExitText = 'Length of a working array is too short. Increase print level.';
   case 6
      ExitText = 'There are false dimensions.';
   case 7
      ExitText = 'The search direction is close to zero. Still infeasible';
   otherwise
      ExitText = ['Unknown Error ' num2str(Inform)];
end

if Inform > 10
      ExitText = ['QP sub problem error ' num2str(Inform-10)];
end
   
Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 030601 hel  Created DFNLP interface.
% 040113 med  Setting up interface.
% 040209 med  ksNLP_x reset
% 040407 ango Also ksNLP_xJ reset, mode parameter removed. Changes to dfgrad.
% 041202 hkh  Remove extra bad call to defblbu
% 041202 hkh  Revise calls to defblbu and StateDef, avoid vector reshuffling
% 041202 hkh  Check if NLLS or approximation problem
% 041202 hkh  SolvType should be 6, not 4