% TOMLAB NLPLSQ Constrained Nonlinear Least Squares Solver
%
% function Result = nlplsqTL(Prob)
%
% INPUT:
%
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% PriLevOpt Print level in solver (0-4)
%
% -----------------------------------------------
% Fields used in Prob.NLPLSQ:
% -----------------------------------------------
%
% ressize   The user must indicate a guess for the approximate size of the
%           objective funtion.  RESSIZE must not be negative.
%
% maxfun    The integer variable defines an upper bound for the number
%           of function calls during the line search.
%
% maxit     Maximum number of outer iterations, where one iteration
%           corresponds to one formulation and solution of the quadratic
%           programming subproblem, or, alternatively, one evaluation of
%           gradients.
%
% maxnm     Stack size for storing merit function values at previous 
%           iterations for non-monotone line search (e.g 10).
%
% tolnm     Relative bound for increase of merit function value, if line
%           search is not successful during the very first step.
% 
% acc       The user has to specify the desired final accuracy (e.g. 
%           1.0e-7).  The termination accuracy should not be smaller than
%           the accuracy by which gradients are computed.
%
% accqp     The tolerance is needed for the QP solver to perform several 
%           tests, for example whether optimality conditions are satisfied
%           or whether a number is considered as zero or not.  If ACCQP is
%           less or equal to zero, then the machine precision is computed 
%           by NLPQLP and subsequently multiplied by 1.0e4
%
% PrintFile Name of NLPLSQ Print file.  Amount and type of printing 
%           determined by PriLevOpt.
%
% ------------------------------------------------------------------------
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
% ExitFlag Exit status.
% ExitText Exit text from solver.
%
% Inform   NLPJOB information parameter.
%
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% Solver   Name of the solver (NLPLSQ).
% SolverAlgorithm  Description of the solver.
%
% NLPJOB.u    Contains the multipliers with respect to the actual iterate
%             stored in x_k.  See user's guide.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc., $Revision: 1.1 $
% Written Jul 27, 2009.    Last modified Jul 27, 2009.

function Result = nlplsqTL(Prob)

Prob.solvType = 6;          % Constrained LS solver, cls
Prob = iniSolve(Prob,3,1,1);% Init globals

Prob = ksDef(Prob);         % Klaus Schittkowski problem definitions

if ~isfield(Prob,'NLPLSQ')  % Make sure NLPLSQ field exists
   Prob.NLPJOB = [];
end

Result=ResultDef(Prob);     % Define the Result struct
Result.Solver='NLPLSQ';     % Solver name
Result.SolverAlgorithm='Dense SQP Maximum-Norm';
BIG=1E12;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

m = m1+m2;
mEQ = Prob.mEQ;

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

% Initial checks on the inputs
xl = bl(1:n); % Prob.x_L
xu = bu(1:n); % Prob.x_U

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( xl,min(xu,x_0(:) ) ); 

[r,c] = nlresid(x_0, Prob);
L = length(r);  % Number of objective functions (or residuals)

% Define default solver options.
NLPLSQ   = DefPar(Prob,'NLPLSQ',[]);
acc       = DefPar(NLPLSQ,'acc',1e-14);
accqp     = DefPar(NLPLSQ,'accqp',1e-14);
ressize   = DefPar(NLPLSQ,'ressize',0.0);
maxfun    = DefPar(NLPLSQ,'maxfun',20);
maxit     = DefPar(NLPLSQ,'maxit',100);
maxnm     = DefPar(NLPLSQ,'maxnm',0);
tolnm     = DefPar(NLPLSQ,'tolnm',0.0);
options   = [acc, accqp, ressize, maxfun, maxit, maxnm, tolnm];
PriLevOpt = DefPar(Prob,'PriLevOpt',0);
PrintFile = DefPar(NLPLSQ,'PrintFile','nlplsq.txt');

[x, Inform, Iter, cLamda, Obj, f_eval, g_eval] = nlplsq ( ...
    L, m, mEQ, x_0, xl, xu, options, PriLevOpt, PrintFile, Prob);
 
Result.f_k  = Obj;
Result.c_k  = [];
Result.cJac = [];
Result.x_k  = x;

optParam = Prob.optParam;

Result = StateDef(Result, x(1:n), [], Result.c_k, ...
   optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);

Result.NLPLSQ.u = cLamda;

ExitFlag = 99;
switch(Inform)
    case -2        
        ExitText = 'Unknown Error: -2';
    case -1
        ExitText = 'Unknown Error: -1';
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
        ExitText = 'The line search could not be terminated successfully.';
    case 5
        ExitText = 'Length of a working array is too short. Increase print level.';
    case 6
        ExitText = 'There are false dimensions.';
    case 7
        ExitText = 'The search direction is close to zero. Still infeasible';
    case 8
        ExitText = 'The starting point violates a lower or upper bound';
    case 9
        ExitText = 'Wrong input parameter.';
    case 10
        ExitText = 'Internal inconsistency of the quadratic subproblem.';
    case {Inform > 100}
        ExitText = ['QP sub problem error. Constraint #' num2str(Inform-100) ' inconsistent.'];
    otherwise
        ExitText = 'Unknown Error';
end

Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 090727 bjo Created NLPLSQ interface.
