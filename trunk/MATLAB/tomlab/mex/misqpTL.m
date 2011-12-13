% TOMLAB MISQP Solver
%
% function Result = misqpTL(Prob)
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
% Prob.MIP      Structure with fields defining the integer properties of
%               the problem. The following fields are used:
%
%   IntVars:
%               If empty, all variables are assumed non-integer
%               If islogical(IntVars) (=all elements are 0/1), then
%               1 = integer variable, 0 = continuous variable.
%               If any element >1, IntVars is the indices for integer variables
%
% -----------------------------------------------
% Fields used in Prob.MISQP:
% -----------------------------------------------
%
% maxit     Maximum number of outer iterations, where one iteration
%           corresponds to one formulation and solution of the quadratic
%           programming subproblem, or, alternatively, one evaluation of
%           gradients.
%
% maxpen    Maximum number of successive increments of the penalty parameter
%           without success (e.g. 50).
%
% maxund    Maximum number of successive iterations without improvements
%           of the iterate X (e.g. 10).
%
% maxnde    Maximum number of branch-and-bound steps for solving MIQP.
%
% nonmon    Maximum number of successive iterations, which are to be
%           considered for the non-monotone trust region algorithm.
%
% brrule    Branching rule.
%           1  Maximal fractional branching (default).
%           2  Minimal fractional branching.
%
% nsrule    Node selection rule.
%           1  Best of all.
%           2  Best of two.
%           3  Depth first (default).
%
% resopt    If set to a nonzero value and if integer variables exist, an additional
%           restart will be performed to check whether an improvement
%           is possible (recommended for functions with curved narrow valleys).
%
% scale     MISQP internally scales continuous variables.
%
% bmod      MISQP modifies the Hessian approximation in order to get more
%           accurate search directions. Calculation time is increased in 
%           case of integer variables.
%
% posort    Transform problem to positive orthant (MIQL).
%
% maxsws    Maximum number of successive warmstarts (default 100).
%
% impbnd    Set nonzero to calculate improved bounds when Best-of-All 
%           node selection strategy is used (default 0).
%
% dfdir     Set nonzero to select direction for Depth-First according to 
%           value of Lagrange function (default 0).
%
% cholmd    Cholesky decomposition mode. 
%           0  Calculate Cholesky decomposition once and reuse it.
%           1  Calculate new Cholesky decomposition whenever warmstart is
%              not active (default).
%
% cutcon    Control the cutting process
%           0  No cuts (default).
%           1  Use disjunctive cuts only.
%           2  Complemented mixed integer rounding (CMIR) cuts only.
%
% maxdcr    Maximal number of rounds for disjunctive cuts, if cutcon > 0.
%
% maxccr    Maximal number of rounds for CMIR cuts.
%
% heurmd    Primal heuristic mode
%           0  No primal heuristic (default).
%           1  Nearest integer.
%           2  Feasibility pump.
%
% acc       The user has to specify the desired final accuracy (e.g. 1.0e-7).
%           The termination accuracy should not be smaller than the accuracy
%           by which gradients are computed. If acc<=0, the machine precision is
%           computed by MISQP and subsequently multiplied by 1.0E4.
%
% accqp     The tolerance is needed for the QP solver to perform several
%           tests, for example whether optimality conditions are satisfied
%           or whether a number is considered as zero or not. If ACCQP is
%           less or equal to zero, then the machine precision is computed
%           by MISQP and subsequently multiplied by 1.0E4.
%
% peninc    Factor for increasing a penalty parameter (SIGMA?), must be greater than
%           one (default 11).
%
% delinc    Factor for increasing the internal descent parameter DELTA.
%           Must be greater than one (default 11).
%
% sigma     Penalty parameter (default 10.0)
%
% delta     Scaling constant (default 10.0)
%
% itrc      Initial trust region radius for continuous variables (default 10.0)
%
% itri      Initial trust region radius for integer variables (default 10.0)
%
% PrintFile Name of MISQP Print file. Amount and type of printing determined
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
% Inform   MISQP information parameter.
%
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% Solver   Name of the solver (MISQP).
% SolverAlgorithm  Description of the solver.
%
% MISQP.u     On return, the first M locations contain the multipliers of the
%             M nonlinear constraints, the subsequent N locations the multipliers
%             subject to the lower bounds, and the final N locations the multipliers
%             subject to the upper bounds. Note that MISQP translates TOMLAB's
%             linear/nonlinear constraints to an internal format, upon which the
%             contents of MISQP.u is based.
%
% MISQP.c     On return, c contains the last computed approximation of the Hessian
%             matrix of the Lagrangian function.
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2007 by Tomlab Optimization Inc., $Release: 5.8.0 $
% Written Sep 1, 2006.    Last modified Feb 23, 2007.

function Result = misqpTL(Prob)

if nargin < 1, error('misqpTL needs the Prob structure as input');end

Prob.solvType = checkType('minlp');          % NLP (CON) solver
Prob = iniSolve(Prob,3,1,1);% Init globals

Prob = ksDef(Prob);         % Klaus Schittkowski problem definitions

if ~isfield(Prob,'MISQP')   % Make sure MISQP field exists
   Prob.MISQP = [];
end

Result=ResultDef(Prob);     % Define the Result struct
Result.Solver='MISQP';      % Solver name
Result.SolverAlgorithm='Dense Mixed Integer SQP';

PriLev=Prob.PriLevOpt;      % Printing level in solver

%
% Define lower and upper bound arrays for MISQP
%
% Inf are changed to BIG (=1E12), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints
%

BIG=1E12;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

% Initial checks on the inputs

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

xl = bl(1:n); % Prob.x_L
xu = bu(1:n); % Prob.x_U

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( xl,min(xu,x_0(:) ) ); 

Result.x_0 = x_0;

mTot = Prob.m;
me  = Prob.mEQ;
n   = size(x_0,1);

[f, g ] = nlfunc(x_0, Prob);
[df,dg] = nlgrad(x_0, Prob);
Result.f_0 = f;

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers, stored in double
IV = zeros(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) || all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('misqp: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
% IntVars = find(IV);
% Detect binary variables (if any integers exist at all)
if any(IV)
   ix = (IV>0 & (bl(1:n)==0.0 & bu(1:n)==1.0));
   IV(ix) = 2;
end

% Find indices of linear equalities and inequalities.
if m1 > 0
  linear = [find(bl(n+1:n+m1)==bu(n+1:n+m1)); ...
            find(bl(n+1:n+m1)~=bu(n+1:n+m1))];
else
  linear = [];
end

PrintFile = DefPar(Prob.MISQP,'PrintFile','misqp.txt');

% Define default solver options.
MISQP = DefPar(Prob,'MISQP',[]);
maxit  = DefPar( MISQP, 'maxit',  []);
maxpen = DefPar( MISQP, 'maxpen', 50);
maxund = DefPar( MISQP, 'maxund', 10);
nonmon = DefPar( MISQP, 'nonmon', 2);
brrule = DefPar( MISQP, 'brrule', 1);
nsrule = DefPar( MISQP, 'nsrule', 3);
% minimum limit on nodes for recomputing Lagrange multipliers
maxnde = DefPar( MISQP, 'maxnde', (2*n+m1+m2+4)^2); 
% minimum limit on nodes for guaranteed generation of disjunctive cuts
% maxnde = DefPar( MISQP, 'maxnde', (2*n+2*(m1+m2)+6)^2); 
resopt = DefPar( MISQP, 'resopt', 1);
scale  = DefPar( MISQP, 'scale' , 0);
bmod   = DefPar( MISQP, 'bmod'  , 0);
posort = DefPar( MISQP, 'posort', 0);
maxsws = DefPar( MISQP, 'maxsws', 100);
impbnd = DefPar( MISQP, 'impbnd', 0);
dfdir  = DefPar( MISQP, 'dfdir' , 0);
cholmd = DefPar( MISQP, 'cholmd', 1);
cutcon = DefPar( MISQP, 'cutcon', 0);
maxdcr = DefPar( MISQP, 'maxdcr', 1);
maxccr = DefPar( MISQP, 'maxccr', 1);
heurmd = DefPar( MISQP, 'heurmd', 0);

% For selected fields: if not set by the user, get values from optParam: 
optParam = DefPar(Prob,'optParam');
if isempty(maxit)
   maxit = DefPar(optParam,'MaxIter',200);
end

acc    = DefPar( MISQP, 'acc'     ,1e-6);
accqp  = DefPar( MISQP, 'accqp'   ,1e-6);
peninc = DefPar( MISQP, 'peninc'  ,11  );
delinc = DefPar( MISQP, 'delinc'  ,11  );
sigma  = DefPar( MISQP, 'sigma'   ,10.0);
delta  = DefPar( MISQP, 'delta'   ,0.01);
itrc   = DefPar( MISQP, 'itrc'    ,10.0);
itri   = DefPar( MISQP, 'itri'    ,10.0);

% Could be included as an option in Prob.MISQP
% ideriv is a vector of length n specifying if the integer variables have 
% derivatives provided by the user ( ideriv(i) = 1 ), or if they should 
% be estimated internally ( ideriv(i) = 0 ), where i is the index of an 
% integer/binary variable. 
% The positions corresponding to the continuous variables are ignored.
ideriv = ones(n,1);

options = [ maxit , maxpen , maxund , nonmon , brrule , nsrule , ...
      maxnde, resopt, scale , bmod , posort , maxsws , impbnd , dfdir , ...
      cholmd , cutcon , maxdcr , maxccr , heurmd , acc , accqp , peninc , ...
      delinc , sigma , delta , itrc  , itri  ];
   
[x,Inform,f_k,g_k,df_k,dg_k,u,c] =  ...
    misqp(mTot, n, me, f, g, df, dg, x_0, xl, xu, IV, linear, ideriv, ...
     PriLev, PrintFile, options, Prob);

Result.f_k  = f_k;
Result.g_k  = df_k;

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
   optParam.xTol, optParam.bTol, optParam.cTol, bl, bu,1);

Result.MISQP.u = u;
Result.MISQP.c = c;

ExitFlag = 99;
switch(Inform)
    case 0
        ExitFlag = 0;
        ExitText = 'The optimality conditions are satisfied.';
    case 1 % iterlimit
        ExitFlag = 1;
        ExitText = 'The algorithm has been stopped after MAXIT iterations.';
    case 2
        ExitText = 'More than MAXUND iterations for fixing trust region';
    case 3
        ExitText = 'More than MAXPEN updates of penalty parameter';
    case 4
        ExitText = 'Termination at infeasible iterate';
    case 5
        ExitText = 'Termination with zero trust region for integer variables';
    case 6
        ExitText = 'Length of a working array is too short';
    case 7
        ExitText = 'There are false dimensions.';
    case 8 
        ExitText = 'Lower and upper bound of an integer variable inconsistent';
    case 9
        ExitText = 'Linear constraints are inconsistent';
    case 11
        ExitText = 'QL could not solve the QP after 40*(N+M) iterations';
    case 12
        ExitText = 'The termination accuracy is insufficient';
    case 13
        ExitText = 'QL terminated due to internal inconsistency, division by zero';
    case 14
        ExitText = 'Numerical instabilities in QL';
   otherwise
        ExitText = ['Unknown Error ' num2str(Inform)];
end

if Inform > 90
   ExitText = ['QP sub problem error ' num2str(Inform-100)];
   ExitFlag = 10;
end

Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 060919 ango Created MISQP interface.
% 070222 hkh  Revise IntVars handling, use new format
% 070223 hkh  Must not use logical IV, because value of IV is 0,1 or 2
