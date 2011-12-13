% TOMLAB SNOPT8 NLP Solver
%
% function Result = snopt8TL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format.
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% c_L, c_U  Bounds on nonlinear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% PriLevOpt Print level.
% WarmStart If true, use warm start, otherwise cold start. Use WarmDefSOL
%           to set the proper parameters.
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution and slacks from previous run.
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%                 2=superbasic (between bounds), 3=basic (between bounds)
% nS        Number of superbasics from previous run.
% hElastic  Defines which variables are elastic in elastic mode. hElastic(j):
%           0 = variable j is non-elastic and cannot be infeasible.
%           1 = variable j can violate its lower bound.
%           2 = variable j can violate its upper bound.
%           3 = variable j can violate either its lower or upper bound.
% moremem   Add more memory if SNOPT stops with not enough storage message.
%           1E6 is 10MB of memory. Default 0.
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined
%           by SPECS parameters or optPar parameters.
%           Output is written on file snoptpri.txt, if not given.
% SummFile  Name of SOL Summary File.
%           Output is written on file snoptsum.txt, if not given.
% To make snopt to not open and not write anything to file:
%    Set SpecsFile and PrintFile empty.
%    Set optPar(1) = 0 AND optPar(3) = 0.
% optPar    Elements > -999 take effect.
%
% -----------------------------------------------
% How optPar is used for setting SPECS parameters:
% -----------------------------------------------
%
% optPar  Structure with optimization parameters.
%
% SNOPT keywords in optPar(#):
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
% --- The SQP Method I - Printing
% 1.  MAJOR PRINT LEVEL             0        1         111111
%
% --- QP subproblems I - Printing
% 2.  MINOR PRINT LEVEL             0        1         10      0, 1 or 10
%
% --- Frequencies I
% 5.  PRINT FREQUENCY               0        100
% 6.  SUMMARY FREQUENCY             0        100
%
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS OPTIONS LISTING      0        0         1       1 = True
%     Also called SUPPRESS PARAMETERS
%
% --- The SQP Method II - Convergence Tolerances
% 9.  MAJOR FEASIBILITY TOLERANCE   >0       1E-6
%
% --- Nonlinear constraints I
% 10. MAJOR OPTIMALITY TOLERANCE    >0       max(2E-6,(10eps_R)^0.5)
%                                            eps_R == optPar(41)
%     Default relative function precision eps_R gives (10*eps_R)^0.5=1.73E-6;
%
% --- QP subproblems II             Convergence Tolerances
% 11. MINOR FEASIBILITY TOLERANCE   >0       1E-6
% 12. MINOR OPTIMALITY TOLERANCE    >0       1E-6
%
% --- Derivative checking
% 13. VERIFY LEVEL                  -1       -1        3       {-1,0,1,2,3}
% 14. START OBJECTIVE CHECK AT COL  0        1         nnObj
% 15. STOP OBJECTIVE CHECK AT COL   0        nnObj     nnObj
% 16. START CONSTRAINT CHECK AT COL 0        1         nnJac
% 17. STOP CONSTRAINT CHECK AT COL  0        nnJac     nnJac
%
% --- QP subproblems III
% 18. SCALE OPTION                  0        0 or 2    2       2 if LP,0 if NLP
%     Option 1 changed to QP Scaling 0 by SNOPT if no linear constraints
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
%     (The value of SCALE PRINT is not shown in the SNOPT print file)
% 21. CRASH TOLERANCE               0        0.1       <1
%
% --- The SQP Method III
% 22. LINESEARCH TOLERANCE          >0       0.9       <1
%
% --- LU I
% 23. LU FACTORIZATION TOLERANCE    1        100 or 3.99       100 if LP
% 24. LU UPDATE TOLERANCE           1        10  or 3.99       10  if LP
% 25  LU SWAP TOLERANCE             >0       1.22E-4           eps^(1/4)
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
%
% --- QP subproblems IV
% 27. PIVOT TOLERANCE               >0       2.22E-15          10*eps
% 28. CRASH OPTION                  0        3         3       {0,1,2,3}
% 29. ELASTIC WEIGHT                0        10000.0
% 30. ITERATION LIMIT               0        10000             or 20m, if more
%     Maximal sum of minor iterations
% 31. PARTIAL PRICE                 0        10 or 1           10 for LP
%
% --- The SQP Method IV
% 32. MAXIMIZE                      0        0         1       0=min,1=maximize
% 33. FEASIBLE POINT                0        0         1       1=feasible pnt
%
% --- Nonlinear constraints I
% 34. VIOLATION LIMIT               >0       1E6
%
% --- The SQP Method V
% 35. MAJOR ITERATIONS LIMIT        >0       max(1000,3*max(n,m))
% 36. MINOR ITERATIONS LIMIT        >0       500
%     Maximal number of minor iterations, in QP/simplex subproblem
% 37. MAJOR STEP LIMIT              >0       2
%
% --- Hessian Approximation I
% 38. HESSIAN FREQUENCY             >0       99999999
%
% --- The SQP Method VI
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     Is set by snoptTL dependent on Prob.ConsDiff, Prob.NumDiff
% 40. DERIVATIVE LINESEARCH         0        1         1       0=NONDERIVATIVE
%     -  0 is quadratic - gives quadratic, without gradient values
%     -  1 is cubic - gives cubic, always using gradient values
%     Default: 0 if numerical derivatives, otherwise 1
% 41. FUNCTION PRECISION            >0       3.0E-13           eps^0.8=eps_R
% 42. DIFFERENCE INTERVAL           >0       5.48E-7           eps^0.4
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.70E-5           eps^{0.8/3}
% 44. PROXIMAL POINT METHOD         1        1         2       {1,2}
%     Minimize the 1-norm (or 2-norm) of |x-x_0| to find an initial point
%     that is feasible subject to simple bounds and linear constraints
% 45. UNBOUNDED STEP SIZE           >0       1E20
% 46. UNBOUNDED OBJECTIVE           >0       1E15
%
% --- Hessian Approximation II
% 47. HESSIAN FULL MEMORY           0        1         1       =1 if nnL <= 75
% or  HESSIAN LIMITED MEMORY                                   =0 if nnL >  75
%
% --- The SQP Method VII
% 48. SUPERBASICS LIMIT             >0       1+nnL     n       =1 if LP problem
%     SUPERBASICS LIMIT >= REDUCED HESSIAN always, see optPar(71)
%
% --- Hessian Approximation III
% 49. HESSIAN UPDATES               >0       20
%     Maximum number of QN (Quasi-Newton) updates
%     If HESSIAN FULL MEMORY, default is 99999999, otherwise 20
% 50. HESSIAN FLUSH                 >0       99999999
%
% --- Frequencies II
% 51. CHECK FREQUENCY               >0       60
% 52. EXPAND FREQUENCY              >0       10000
% 53. FACTORIZATION FREQUENCY       >0       50
%
% --- LU II
% 63. LU PARTIAL  PIVOTING          0        0         3  0=partial
%     or LU COMPLETE PIVOTING                             1=complete
%     or LU ROOK PIVOTING                                 2=rook
%     or LU DIAGONAL PIVOTING                             3=diagonal
%
% --- The SQP Method VIII
% 64. PENALTY PARAMETER           >=0       0.0
%     Initial penalty parameter
%
% --- QP subproblems V
% 65. NEW SUPERBASICS              >0       99            also MINOR SUPERBASICS
%     Maximal number of new superbasics per major iteration
%
% 66. QPSOLVER CHOLESKY             0        0         2  0=Cholesky
%     or QPSOLVER   CG                                    1=CG
%     or QPSOLVER   QN                                    2=Quasi-Newton CG
%
% --- Conjugate-Gradient QP solver
% 67. CG TOLERANCE                 >0      1E-2
% 68. CG ITERATIONS                >0      100            Max number of CG iters
% 69. HESSIAN PRECONDITIONING       0        0         1  QN preconditioned CG
%     also called
%     QG      PRECONDITIONING       Default 1 if QPSOLVER QN
% 70  SUBSPACE                      0      0.1         1  Subspace tolerance
%     Quasi-Newton QP rg tolerance
%
% --- The SQP Method IX
% 71. HESSIAN DIMENSION            >0  min(2000,nnL+1) n  =1 if LP problem
%     or REDUCED HESSIAN           Number of columns in Reduced Hessian
%
% -----------------------------------------------------------------------
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial  solution vector.
% g_k      Gradient of the function.
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from snopt.m (similar to TOMLAB).
% Inform   SNOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (snopt).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field:
% xs       Solution and slack variables.
% hs       State for variables and slacks in xs.
% nS       #   of superbasics.
% nInf     #   of infeasibilities.
% sInf     Sum of infeasibilities.
%
% -----------------------------------------------------------------------
%
% For a problem description, see conAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written July 4, 2000.   Last modified March 20, 2008.

% -----------------------------------------------
% Experimental parameters, only available using SPECS file:
% -----------------------------------------------
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
% --- QP subproblems
% xx. ELASTIC OBJECTIVE             1        2         2       2=Use elastic
% xx. ELASTIC WEIGHTMAX             0        100.0             Elastic weightmax
% xx. ELASTIC MODE                  0        100.0             >0 use elastic
%
% --- LU
%     Lmax1  = Maximum multiplier allowed in  L  during refactorization.
%     Lmax2  = Maximum multiplier allowed during updates.
% xx. LU DEFAULTS TPP, tolDpp      >0        100 or 3.99       100 for LP
%     Default Lmax1 for TPP, Threshold Partial Pivoting
% xx. LU DEFAULTS TCP, tolDcp      >0        10 or 20          10  for LP
%     Default Lmax1 for TCP, Threshold Complete Pivoting
% xx. LU DEFAULTS UPDATES, tolDup  >0        10 or 3.99        10  for LP
%     Default Lmax2
% xx. LU DENSITY TOLERANCE          >0       0.5         1
%
% --- The SQP Method
% xx. DERIVATIVE OPTION             0        1         2       {0,1,2}
%     Influences if to use derivatives in the line search
%     =0 Assumes derivatives unknown, =1 Assumes known
%     snJac sets DERIVATIVE OPTION = 0 hard coded
%     DERIVATIVE LINESEARCH / NONDERIVATIVE LINESEARCH overrides this parameter
%     NOTE: DERIVATIVE OPTION, not DERIVATIVE LEVEL is written in Print File

% HKH:
% MinorIter has been changed from max(1000,5*max(m,n)) to 500
% (but still used in SQOPT)

function Result = snopt8TL(Prob)

if nargin < 1, error('snopt8TL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 3; % NLP (CON) solver

% ANGO - Need to make this more intelligent. If 2nd derivs are available,
% we should do iniSolve(Prob,3,2,2)

if isempty(Prob.NumDiff)
   ObjDers = 1;
elseif Prob.NumDiff==6
   ObjDers = 0;
else
   ObjDers = 1;
end
if Prob.mNonLin > 0
   if isempty(Prob.ConsDiff)
      ConsDers = 1;
   elseif Prob.ConsDiff==6
      ConsDers = 0;
   else
      ConsDers = 1;
   end
else
   ConsDers = 0;
end

Prob = iniSolve(Prob,3,ObjDers,ConsDers);

Result=ResultDef(Prob);
Result.Solver='snopt';
Result.SolverAlgorithm='SNOPT 8.0-0(0) NLP code';

PriLev=Prob.PriLevOpt;

% Define lower and upper bound arrays for SNOPT
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

BIG=1E20;

% Request bl/bu in order (x,nonlinear constraints, linear constraints)
[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 0);

% Initial checks on the inputs
ObjAdd  = 0;

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isempty(moremem), moremem = 0; end

m       = m1+m2;
nnObj = n; % number of nonlinear variables in objective

% Check if any linear part
probType = Prob.probType;
%if ~any(probType==[2 7 8]) & isfield(Prob.QP,'c')
if isfield(Prob.QP,'c')
   if probType == 2
      c  = [];
      m3 = 0;
      % Does not speed up much to avoid linear part in callback
      %c  = Prob.QP.c;
      %m3 = length(c) > 0;
      %if m3
      %   Prob.FUNCS.f='qp_fX';
      %   Prob.FUNCS.g='qp_gX';
      %end
   elseif strcmpi(Prob.FUNCS.f,'lp_f')
      c  = Prob.QP.c;
      m3 = length(c) > 0;
      nnObj = 0;
   else
      c  = [];
      m3 = 0;
   end
else
   c  = [];
   m3 = 0;
end

% Check on Warm start, then set hs,xs,nS
if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs;
   hs   = Prob.SOL.hs;
   %xs   = Prob.SOL.xs(1:nb);
   %hs   = Prob.SOL.hs(1:nb);
   nS   = Prob.SOL.nS;
   x_0  = xs(1:n);
else
   % Initial values
   x_0  = Prob.x_0(:);
   Warm = 0;
   nS   = Prob.SOL.nS;
   hs   = Prob.SOL.hs;
   if length(x_0) < n, x_0=zeros(n,1); end
   % Safe guard x_0 
   x_0 = max(bl(1:n),min(bu(1:n),x_0(1:n)));
   xs   = x_0;
end

if m2 > 0
   % Determine the sparse problem structure
   % number of nonlinear constraints   (nnCon)
   nnCon=m2;
   if isempty(Prob.ConsPattern)
      gJac=ones(m2,n);
      nnJac=n;
   else
      gJac=Prob.ConsPattern;
      [ix,iy]=find(gJac);
      %[i,j]=find(Prob.ConsPattern);

      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(Prob.ConsPattern),ix,iy);

      % Number of nonlinear Jacobian variables (nnJac)
      nnJac=size(gJac,2);
   end
   % Must call function, in case global values are needed for constraints
   % Because snopt calls constraints first, then functions
   Result.f_0= nlp_f(x_0, Prob);
   % Sometimes maybe also a call to the gradient might be needed.
   % Do not make such a call now, it is highly unlikely.
else
   nnCon      = 0;
   nnJac      = 0;
   gJac       = [];
   Result.f_0= nlp_f(x_0, Prob);
   %Result.f_0 = c'*x_0;
end

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end 
end 

% Set up the constraint matrix A 
if isempty(gJac)
   A=sparse(Prob.A);
else
   if isempty(Prob.A)
      A = sparse(gJac);
   else
      if nnJac < n
         A = sparse([[gJac,zeros(m2,n-nnJac)];Prob.A]);
      else
         A = sparse([gJac;Prob.A]);
      end
   end
end
if m==0
   % Construct one linear constraint from some bounded variable
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);
   bl(n+1,1)=max(bl(i),-BIG);
   bu(n+1,1)=min(bu(i),BIG);
   m=1;
   A = sparse(zeros(1,n));
   A(1,i)=1;
end

if m3 > 0
   A = [A;c'];
   iObj=m+1;  % The constraint row with linear obj term c is put last
   % Add the bounds for the objective in bl and bu
   bl=[bl;-BIG];
   bu=[bu; BIG];
else
   iObj=0;
end

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('snopt',3,...
%         nnObj, nnJac, size(A,1), Prob);

optPar    = Prob.SOL.optPar(:)';
optPar(length(optPar)+1:Prob.SOL.optParN)=-999;

% Avoid gradient evaluations in line search if numerical derivatives
if optPar(40) < 0 & (Prob.NumDiff > 0 | (Prob.ConsDiff > 0 & m2 > 0) | ...
   isempty(Prob.FUNCS.g) | (isempty(Prob.FUNCS.dc) & m2 > 0) | Prob.CheckNaN)
   % Change default to non-Derivative line search if numerical derivatives
   optPar(40) = 0;
end
% SNOPT has changed the logic for SUPERBASICS
% Now it is the Reduced Hessian that might be too small, for problems > 2000
%if optPar(48) < 0 & (n-m1-m2) > 450
%   % Increase number of superbasics
%   if n > 5000
%      optPar(48) = max(500,n-m1-m2+200);
%   else
%      optPar(48) = max(500,n+1);
%   end
%end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

% Lagrange multipliers, known info
pi(1:m+m3)=0;

if Prob.NumDiff == 6
   if Prob.ConsDiff == 6
      optPar(39)=0;
   else
      optPar(39)=2;
   end
elseif Prob.ConsDiff == 6
   optPar(39)=1;
else
   optPar(39)=3;
end

% optPar(39)=4 activates SNOPT 8 asking for exact Hessian*vector products
optPar(39) = 4;
% Avoid derivative check
if optPar(13) < 0
   optPar(13)=-1;
end

[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     snopt8( A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            Prob.SOL.SpecsFile, Prob.SOL.PrintFile, Prob.SOL.SummFile, ...
            PriLev, ObjAdd, moremem, Prob.Name );

Result.f_k=Obj;
Result.g_k=gObj;
Result.c_k=fCon;
Result.x_k=xs(1:n);
Result.x_0=x_0;
% Result.v_k=[rc(1:n);pi(1:m)];
% Multipliers for bounds-linear-nonlinear
Result.v_k=[rc(1:n);pi(m2+1:m2+m1);pi(1:m2)];

if ~isempty(gCon)
   [ix,iy] = find(gJac);
   % TOMLAB now has the format one constraint / row, same as SOL solvers
   Result.cJac=sparse(ix,iy,gCon,size(gJac,1),size(gJac,2));
end

% Saved for a warm start
%if m3 > 0
   Result.SOL.xs=xs;
   Result.SOL.hs=hs;
%else
%   % The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
%   Result.SOL.xs=[xs;-Obj];
%   Result.SOL.hs=[hs;3];
%end
Result.SOL.nS=nS;
Result.SOL.nInf=nInf;
Result.SOL.sInf=sInf;
Result.SOL.optPar=optPar;

Result.FuncEv    = [];%iwCount(3);
Result.GradEv    = [];%sum(iwCount(4:6));
Result.ConstrEv  = [];%iwCount(7);
Result.Iter      = iwCount(2);
Result.MinorIter = iwCount(1);

optParam = Prob.optParam;
if m1 > 0
   Result.Ax = A(m2+1:m2+m1,:)*xs(1:n);
   Result = StateDef(Result, xs(1:n), Result.Ax, fCon, ...
                  optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 0);
else
   Result = StateDef(Result, xs(1:n), [], fCon, ...
                  optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 0);
end

% ExitText
switch(Inform)
   case 0  , Text = 'Finished successfully';
   case 1  , Text = 'Optimality conditions satisfied';
   case 2  , Text = 'Feasible point found';
   case 3  , Text = 'Requested accuracy could not be achieved';
      
   case 10 , Text = 'The problem appears to be infeasible';
   case 11 , Text = 'Infeasible linear constraints';
   case 12 , Text = 'Infeasible linear equalities';
   case 13 , Text = 'Nonlinear infeasibilities minimized';
   case 14 , Text = 'Infeasibilities minimized';
      
   case 20 , Text = 'The problem appears to be unbounded';
   case 21 , Text = 'Unbounded objective';
   case 22 , Text = 'Constraint violation limit reached';
      
   case 30 , Text = 'Resource limit error';
   case 31 , Text = 'Iteration limit reached';
   case 32 , Text = 'Major iteration limit reached';
   case 33 , Text = 'Superbasics limit is too small';
      
   case 40 , Text = 'Terminated after numerical difficulties';
   case 41 , Text = 'Current point cannot be improved';
   case 42 , Text = 'Singular basis';
   case 43 , Text = 'Cannot satisfy the general constraints';
   case 44 , Text = 'Ill-conditioned null-space basis';
      
   case 50 , Text = 'Error in the user-supplied functions';
   case 51 , Text = 'Incorrect objective derivatives';
   case 52 , Text = 'Incorrect constraint derivatives';
      
   case 60 , Text = 'Undefined user-supplied functions';
   case 61 , Text = 'Undefined function at the first feasible point';
   case 62 , Text = 'Undefined function at the initial point';
   case 63 , Text = 'Unable to proceed into undefined region';
      
   case 70 , Text = 'User requested termination';
   case 72 , Text = 'Terminated during constraint evaluation';
   case 73 , Text = 'Terminated during objective evaluation';
   case 74 , Text = 'Terminated from monitor routine';
      
   case 80 , Text = 'Insufficient storage allocated';
   case 81 , Text = 'Work arrays must have at least 500 elements';
   case 82 , Text = 'Not enough character storage';
   case 83 , Text = 'Not enough integer storage';
   case 84 , Text = 'Not enough real storage';
      
   case 90 , Text = 'Input arguments out of range';
   case 91 , Text = 'Invalid input argument';
   case 92 , Text = 'Basis file dimensions do not match this problem';
   case 140, Text = 'System error';
   case 141, Text = 'Wrong number of basic variables';
   case 142, Text = 'Error in basis package';
   otherwise
      Text = 'NOTE: UNKNOWN SNOPT Inform value.';
end

% ExitFlag
if(Inform<10)     % Solution found, not necessarily to requirements
  ExitFlag = 0;
elseif(Inform<20) % Infeasible
  ExitFlag = 4;
elseif(Inform<30) % Unboundedness
  ExitFlag = 2;
elseif(Inform<40) % Resource limit (iterations etc)
  ExitFlag = 1;
else              % Input error. Not very accurate.
  ExitFlag = 10;
end

Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;
Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nSNOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('SNOPT: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Major       iterations%7d. ',iwCount(2));
   fprintf('Total minor iterations%7d. ',iwCount(1));
   fprintf('\n');
   
   fprintf('fObj and gObj evaluations%7d %d %d %d\n',iwCount(3:6));
   if m2 > 0
      fprintf('fCon and gCon evaluations%7d %d %d %d\n',iwCount(7:10));
   end
   fprintf('nInf (# of infeasible constraints) %7d. ',nInf);
   fprintf('nS (# of superbasics) %7d. ',nS);
   fprintf('sInf (Sum of infeasibilities outside bounds) %14.7e\n',sInf);

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(xs(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector hs for x and constraints = \n');
      xprinti(hs(1:min(length(hs),MAX_x)),'hs: ');
   end
   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (pi) = \n');
      xprinte(pi(1:min(length(pi),MAX_c)),'pi:');
      fprintf('Reduced costs rc: Last %d elements should be v_k\n',m);
      xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
   end
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000704 hkh New snoptMex interface.
% 000916 hkh Return string ExitText with interpretation of Inform flag 
% 010901 hkh Give snopt correct number of nonlinear Jacobian variables 
%            using size of Prob.ConsPattern
% 010903 hkh Add complete pivoting option. Allow cold start change of hs,ns
% 011212 hkh Improve comments on iterations
% 011212 hkh Use Prob.Name instead of creating new variable ProbName
% 020210 hkh Send linear index for sparse gJac subscripts in Prob.ConsIdx
% 020304 hkh Set optPar(39) DERLVL dependent on Prob.ConsDiff, Prob.NumDiff
% 020417 hkh Send optPar back in Result.SOL
% 021105 hkh Make a comment about SCALE OPTION 1 set to 0 when only linear con
% 030801 hkh Correct test for LP problem, avoid f==0 results
% 031208 ango Edit comments (optPar(63)), Result.SolverAlgorithm for SNOPT 6.2-2
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040602 med  Help fix PriLev to PriLevOpt
% 041203 hkh  Removed unnecessary print statements from old debug
% 041203 hkh  Revise calls to defblbu and StateDef, use different Order=0
% 041203 hkh  If internal differentiation, different call to iniSolve 
% 041203 hkh  Added safe guard for x_0, needed
% 041203 hkh  Change order in Result.v_k: bounds, linear, nonlinear
% 050212 med  Default for scaling corrected in help
% 050318 ango Fix ExitFlag
% 050605 hkh  Correction of defaults for optPar 2,10,27,28,29,34,35,36,42
% 050605 hkh  New 6.2 parameter LU SWAP TOLERANCE as optPar(25), skip LU DENSITY
% 050605 hkh  % 25  LU DENSITY TOLERANCE          >0       0.6 (NOT USED in code)                            
% 050605 hkh  New 6.2 parameter PROXIMAL POINT METHOD as optPar(64)
% 050605 hkh  New 6.2 parameter PENALTY PARAMETER as optPar(65)
% 050605 hkh  New 6.2 parameter NEW SUPERBASICS (MINOR SUPERBASICS) as optPar(66)
% 050605 hkh  % 44. FEASIBLE EXIT   0  0  1  1=feasible exit (Obsolete)
% 050606 hkh  Adding comments about parameters only available using SPECS
% 050606 hkh  New 7.1 parameter: 66. QPSOLVER CHOLESKY,CG or QNCG
% 050606 hkh  New 7.1 parameters: 67. CG TOLERANCE, 68. CG ITERATIONS
% 050606 hkh  New 7.1 parameter: HESSIAN (or CG) PRECONDITIONING
% 050606 hkh  New 7.1 parameter: SUBSPACE
% 050606 hkh  New 7.1 parameter: 71. HESSIAN DIMENSION (also REDUCED HESSIAN)
% 050606 hkh  Removed special setting of SUPERBASICS (optPar(48))
% 050614 med  Updated help for optPar vector
% 050614 med  Removed optPar 3 and 4, 54-62
% 050725 med  bl, bu safe-guarded to be a column vector
% 050827 hkh  If CheckNaN=1, non-derivative line search
% 050829 med  FuncEv, GradEv, ConstrEv set to []
% 050927 med  Help updated for QN subsolver option
% 051216 med  optPar(13) default to -1 (no derivative check)
% 060814 med  FUNCS used for callbacks instead
% 061025 med  No callbacks anymore if objective is linear
% 061224 med  moremem help updated
% 080320 hkh  Avoid optPar(40)=0 if Prob.NumDiff < 0 | Prob.ConsDiff < 0
% 080320 hkh  Removed SOLSet call
