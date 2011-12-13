% TOMLAB LGO nonlinear constrained global optimization solver
%
% function Result = lgoTL(Prob)
%
% LGO solves continuous global constrained nonlinear problems of the
% type:
%
%    min f(x)
%
%    subject to
%
%    x_L <= x <= x_U
%
%    g(x) <= 0
%    h(x) == 0
%
% where x,x_L,x_U are n-vectors and [ g(x);h(x) ] is an m-vector of
% linear and/or nonlinear functions.
%
% ---------------------------------------------------------------------
%
% The standard TOMLAB global constrained nonlinear (glc) problem definition is:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, n  variable bounds
%    b_L <= A*x  <= b_U, m1 linear constraints
%    c_L <= c(x) <= c_L, m2 nonlinear constraints
%
% where x,x_L,x_U are n-vectors,
% b_L,b_U are m1-vectors, A is a m1xn matrix,
% c_L,c_L are m2-vectors and c(x) is an m2-vector of nonlinear functions.
%
% LGO does NOT treat integer variables, which is possible to define
% with the glc
%
% lgoTL automatically handles the translation from the TOMLAB
% formulation to the format needed by LGO, and back.
%
% This is done by the separation of constraints with both upper and
% lower bounds.  For maximum performance, the user is advised to
% formulate his/her problems with single-bounded constraints if
% possible. Each double-bounded constraint will be expanded to TWO
% distinct constraints internally.
%
% -----------------------------------------------------------------------
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%  x_L, x_U     Bounds on variables.
%  b_L, b_U     Bounds on linear constraints.
%  A            Linear constraint matrix.
%  c_L, c_U     Bounds on nonlinear constraints.
%               For equality constraints (or fixed variables), set
%               e.g. b_L(k) == b_U(k).
%
%  PriLevOpt    Print level in the solver.
%
%  MaxCPU       Maximum runtime in seconds. This parameter is overridden by
%               Prob.LGO.options.tlimit, if set.
%
%  optParam     Structure with optimization parameters.
%
%               Fields used:
%
%    MaxFunc    Maximum number of function evaluations.
%               (Prob.LGO.options.g_maxfct)
%
%  LGO          Structure with LGO solver specific fields.
%
%    PrintFile  Name of file to print progress information and results to.
%               Default: 'lgo.out'
%
%    SummFile   Name of file to print summary information and results to.
%               Default: 'lgo.sum'
%
%  LGO.options  Structure with special fields for LGO optimization
%               parameters. The following fields are used:
%
%      opmode     Operation mode (Default 3):
%
%            0    Local search (only) from given nominal solution (LS)
%            1    Global branch-and-bound search + local search   (BB+LS)
%            2    Global adaptive random search + local search    (GARS+LS)
%            3    Global multistart random search + local search  (MS+LS)
%
%      irngs      Random generator seed value. (Default 0)
%
%      g_maxfct   Global search termination criterion parameter:
%                 global scope search ends, if the number of merit function
%                 evaluations (approximately) attains g_maxfct.
%                 For G_MAXFCT=0 local search starts from the given nominal
%                 solution. (Default -1)
%
%      max_nosuc  global search termination criterion parameter:
%                 global search phase ends, if the current best solution
%                 did not improve during (at least) the last MAXNOSUC
%                 subsequent merit function evaluations. (Default -1)
%
%      tlimit     Program execution time limit (seconds) (Default 600)
%
%      penmult    Constraint penalty multiplier. Global merit function is defined
%                 as objective + the violated constraints weighted by penmult.
%                 (Default 100)
%
%      acc_tr     Global search termination criterion parameter: global phase
%                 ends, if the overall merit function value found in global
%                 search is less than acc_tr. The merit function is the sum
%                 of the objective function and the violated constraints
%                 (if present). In the global search phase, the latter are added
%                 to the objective, applying unit penalty multiplier factors
%                 (set within LGO). (Default -1.0e10)
%
%      fi_tol     Local search termination criterion parameter:
%                 first local search phase ends, if the merit function
%                 improvement is less than fi_tol. (Default  1.0e-6)
%
%      fct_trg    Large objective function value in local search; partial
%                 stopping criterion in second local search phase.
%                 If unknown, then can be set to 'safe' lower bound, in order
%                 to skip criterion. (Default -1.0e10)
%
%      con_tol    Maximal constraint violation tolerance in local search;
%                 partial stopping criterion in second local search phase.
%                 (Default 1.0e-6)
%
%      kt_tol     Tolerance in satisfying the Kuhn-Tucker local optimality
%                 conditions; stopping criterion in third local search phase.
%                 (Default 1.0e-6)
%
% ----------------------------------------------------------------------------------
%
% OUTPUT:
% Result        Structure with optimization results.
%
%   f_k         Function value at optimum.
%
%   x_k         Solution vector.
%   x_0         Initial solution vector.
%
%   c_k         Nonlinear constraint residuals.
%
%   xState      State of variables. Free == 0; On lower == 1; On upper == 2;
%               Fixed == 3;
%
%   bState      State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   cState      State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   Ax          Values of linear constraints.
%
%   ExitFlag    Exit status.
%   ExitText    Exit text from solver.
%   Inform      LGO information parameter. See Result.ExitText
%
%   FuncEv      Number of function evaluations.
%   ConstrEv    Number of constraint evaluations (same as FuncEv)
%   QP.B        Basis vector in TOMLAB QP standard.
%
%   Solver           Name of the solver (LGO).
%   SolverAlgorithm  Description of the solver.
%
%
%   LGO         Subfield with LGO specific results
%      sstat    Solver status
%      mstat    Model status
%      runtime  Time spent, measured by LGO solver.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Jan 8, 2004.   Last modified Aug 23, 2005.

function Result = lgoTL(Prob)

%#function nlfunc

Prob.solvType = 12;

Prob = iniSolve(Prob,12,0,0);

Result = ResultDef(Prob);
Result.Solver = 'LGO';

BIG = 1E6;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);
sense  = -1;
Prob   = ksDef(Prob,sense);
m      = Prob.m;
ctype  = -1*ones(m,1);
ctype( Prob.AixEQ    ) = 0;
ctype( Prob.cixEQ+m1 ) = 0;

% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)

PriLev = Prob.PriLevOpt;

% Starting point
x0        = DefPar(Prob,'x_0',0.5*(bl(1:n)+bu(1:n)));
% Safe-guarded starting point
x0        = max(bl(1:n),min(x0,bu(1:n)));
LGO       = DefPar(Prob,'LGO',[]);
options   = DefPar(LGO,'options',[]);

irngs     = DefPar(options,'irngs'     ,0       );  
g_maxfct  = DefPar(options,'g_maxfct'  ,-1      );
max_nosuc = DefPar(options,'max_nosuc' ,-1      );
penmult   = DefPar(options,'penmult'   ,100.0   );
fi_tol    = DefPar(options,'fi_tol'    ,1.0e-6  );
con_tol   = DefPar(options,'con_tol'   ,1.0e-6  );
kt_tol    = DefPar(options,'kt_tol'    ,1.0e-6  );

acc_tr    = DefPar(options,'acc_tr'    ,-1.0e10 );
fct_trg   = DefPar(options,'fct_trg'   ,-1.0e10 );

tlimit    = DefPar(options,'tlimit'    ,[]      );

% Prob.LGO.options.tlimit overrides Prob.MaxCPU
if isempty(tlimit)
   if ~isempty(Prob.MaxCPU)
      tlimit = Prob.MaxCPU;
   else
      tlimit = 600;     % LGO default
   end
end

if isinf(tlimit), tlimit = 600; end

ioptions = [ irngs(1) g_maxfct(1) max_nosuc(1) tlimit(1) ];
foptions = [ penmult(1) fi_tol(1) con_tol(1) kt_tol(1)   ];

% First check Prob.Alg for an algorithm flag. Then check if there's
% been something set in Prob.LGO.options.opmode - which takes precedence. 
opmode = DefPar(Prob,'Alg',3);
opmode = DefPar(options,'opmode',opmode);

% Last of all - if illegal value, choose mode=3 (?)
if(opmode<0 | opmode>3), opmode = 3; end

switch(opmode)
 case 0
  Result.SolverAlgorithm='LGO solver (Local mode)';
 case 1
  Result.SolverAlgorithm='LGO solver (Branch & bound mode)';
 case 2
  Result.SolverAlgorithm='LGO solver (Adaptive random search mode)';
 case 3
  Result.SolverAlgorithm='LGO solver (Multistart random search mode)';
end

PrintFile = DefPar(LGO,'PrintFile','lgo.out');
SummFile  = DefPar(LGO,'SummFile','lgo.sum');

[ x_k,f_k,g_k,nfev,sstat,mstat,runtime ] = ...
    lgo( n,m,bl(1:n),bu(1:n),x0,ctype, opmode, ...
	 acc_tr,fct_trg,ioptions, foptions, ...
	 PriLev,PrintFile,SummFile,Prob);

Result.LGO.mstat   = mstat;
Result.LGO.sstat   = sstat;
Result.LGO.runtime = runtime;

c_k    = zeros(m2,1);
cixEQ  = Prob.cixEQ;
cixLow = Prob.cixLow;
cixUpp = Prob.cixUpp;
if ~isempty(cixEQ)
   c_k(cixEQ)  = sense*(g_k(cixEQ) + Prob.c_L(cixEQ));
end
if ~isempty(cixLow)
   c_k(cixLow) = sense*(g_k(cixLow) + Prob.c_L(cixLow));
end
if ~isempty(cixUpp)
   c_k(cixUpp) = sense*(Prob.c_U(cixUpp) - g_k(cixUpp));
end

if m1>0, Result.Ax  = Prob.A*x_k; else Result.Ax = [];  end

Result.x_k = x_k;
Result.f_k = f_k;
Result.c_k = c_k;

xTol = DefPar(Prob.optParam,'xTol',kt_tol);
bTol = DefPar(Prob.optParam,'bTol',con_tol);
cTol = DefPar(Prob.optParam,'cTol',con_tol);

Result = StateDef(Result, x_k, Result.Ax, Result.c_k, xTol, ...
                  bTol, cTol, bl, bu, 1);

% SOLVER STATUS
% 1    Normal completion
% 2    Iteration interrupt (not used; will be implemented)
% 3    Resource interrupt (time limit)
% 4    Terminated by solver (due to set solver parameters or solver limitations)
% 7    Built-in LGO size limitations are exceeded by model

switch(sstat)
 case 1, sText = 'Normal completion';
 case 2, sText = 'Iteration interrupt';
 case 3, sText = 'Time limit exceeded';
 case 4, sText = 'Terminated by solver';
 case 7, sText = 'Size limitation';
 otherwise, sText = 'Unknown solver status';
end
Result.Inform = sstat;

%! MODEL STATUS
% 1 Globally optimal (the solution is generated by using global+local search)
% 2 Locally optimal  (the solution is generated by using only local search)
% 3 Unbounded (in presence of bounds, this should happen only due to modeling errors)
% 4 Infeasible (this can occur also when the solution is not sufficiently accurate)
% 7 Intermediate this value is returned (w/o text string below) upon SST=3;

switch(mstat)
 case 1, mText ='Globally optimal solution found';  ExitFlag = 0;
 case 2, mText ='Locally optimal solution found';   ExitFlag = 0;
 case 3, mText ='Unbounded solution indicated';     ExitFlag = 2;
 case 4, mText ='Infeasible solution found';        ExitFlag = 1;
 case 7, mText ='Intermediate solution returned';   ExitFlag = 0;
 otherwise, mText = 'Unknown model status';         ExitFlag = 10;
end

Result.ExitText = [sText ' : ' mText];
Result.ExitFlag = ExitFlag;

Result=endSolve(Prob,Result);

% MODIFICATION LOG
%
% 040108 ango File created
% 040109 ango Call iniSolve/endSolve
% 040109 hkh  Use m2 for original number of nonlinear constraints
% 040112 ango Fix defaults for acc_tr and fct_trg
% 040122 ango Fix constraint handling, comments. 
% 040122 med  Added help info.
% 040124 hkh  Someone has made a comment on Prob.solvType - error
% 040124 hkh  Default changed for g_maxfct and max_nosuc
% 040728 med  Pragma added for MATLAB Compiler
% 041020 ango Support for MaxCPU added
% 041202 hkh  Revise calls to defblbu and StateDef, avoid vector shuffling
% 041222 med  Safeguard added for x_0
% 050117 med  Options changed to options
% 050503 ango Safeguard for infinite tlimit  
% 050823 ango Change default tlimit from 1e8 to 1e7
% 060131 med  Time limit set to original default (600)