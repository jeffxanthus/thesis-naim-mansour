% multiMin.m
%
% multiMin solves general constrained mixed-integer global nonlinear
% optimization problems. It tries to find all local minima by a multi-start
% method using a suitable nonlinear programming subsolver or nonlinear
% least squares solver, unless the user specifies another solver.
%
% If Parallel Computing Toolbox is installed and Prob.Threads > 1, 
% then multiMin is running 2 major loops in parallel with parfor,
% giving a speedup roughly Prob.Threads-1
%
% By default multiMin uses the TOMLAB license dependent large-scale
% nonlinear programming solver, but the user may change this by setting
% Prob.GO.localSolver. If multimin detects a nonlinear least squares
% problem, unconstrained or constrained, it uses the TOMLAB license
% dependent dense nonlinear least squares solver.
%
% Either multiMin generates m initial points randomly inside the box
% defined by the simple bounds for the problem or uses a matrix of initial
% points given by the user to complete the local searches.
%
% If generating random points and there are linear constraints, multiMin
% checks feasibility with respect to the linear constraints, and for each
% initial trial point tries 100 times to get linear feasible by generating
% new random points before accepting the initial point, even if linear
% infeasible.
%
% If generating random points and Prob.x_0 is non-empty, multiMin uses
% Prob.x_0 as the 1st initial point (even if not feasible), i.e. only m-1
% initial points are random in this case.
%
% If fGoal (Prob.optParam.fGoal) is set, multiMin stops when a local
% optimum with function value less or equal to fGoal is found.
%
% multiMin checks the function value f(x_0) for every initial point x_0,
% and does no local search if f(x_0) > fCut. The user can input fCut
% (Prob.fCut), default 1E10.
%
% multiMin solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% The integer components of every x_0 are rounded to the nearest integer
% value inside simple bounds, and these components are fixed during the
% nonlinear local search.
%
% Note that it is recommended to set the simple bounds x_L and x_U as tight
% as possible, because the random points are sampled anywhere in this box
% The tighter the bounds, the more likely that multiMin will find the global
% minimum.
%
% Usage: See below
%
% Calling syntax:
%
% Result = multiMin(Prob, xInit); or
% Result = tomRun('multiMin', Prob, 1);
%
% INPUT PARAMETERS
%
% xInit     Either
%           1x1 number - 
%           If > 0, the number of random initial points, default min(3000,30*n);
%           If < 0, generate |xInit| points with Latin Hypercube design
%
%           dxm matrix - Every column is an initial point (of length
%           d=Prob.N)
%
%           xInit may also be set as Prob.xInit
%
% Mx0     When doing warm start, multiMin only does local optimization using the "best" Mx0
%         initial x_0 values generated from the input of xInit. "Best" is defined as the
%         x_0 values with largest distance from the X0 and XX points in the previous run.
%
% Prob    Structure, defined as to solve a standard constrained nonlinear
%         problem. The Prob structure is feed to the localSolver.
%         See e.g. conAssign
%
% Fields in the Prob structure used by multiMin itself:
%
%   fCut      If initial f(x_0) > fCut, no local optimization is done
%             Default fCut = 1E10
%
%   WarmStart If true, >0, multiMin assumes that the field Prob.multiMin
%             contains the output from the previous run for the same
%             problem. See the Output field Result.multiMin.
%             Use WarmDefGLOBAL to set the correct fields in the Prob
%             structure. Call Prob = WarmDefGLOBAL('multiMin',Prob,Result);
%
%   RandState If ~WarmStart and isscalar(xInit), 
%             RandState is used as follows:
%             If > 0, rand('state',RandState) is set to initialize the
%             pseudo-random generator
%             if < 0, rand('state',sum(100*clock)) is set to give a new set
%             of random values each run
%             if RandState == 0, rand('state',) is not called
%             Default RandState = -1
%
%   xEqTol    Tolerance to test if new point x_k already defined as
%             optimum: norm(x_k-xOpt(:,i)) <= xEqTol*max(1,norm(x_k))
%             If test fulfilled x_k is assumed to be too close to xOpt(:,i)
%             Default xEqTol = 1E-5
%
%   fEqTol    Tolerance to test if new f(x_k(1)) is same as (x_k(2))
%             Default fEqTol = 1E-7
%
%   x_L       Lower bounds for each element in x. If generating random
%             points, -inf elements of x_L are set to min(-L,xMin,x_U-L)
%             xMin is the minimum of the finite x_L values
%   x_U       Upper bounds for each element in x. If generating random points,
%             inf elements of x_U are set to max(L,xMax,x_L+L)
%             xMax is the maximum of the finite x_U values
%
%             L is 100 for nonlinear least squares, otherwise 1000
%
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%
%   PriLevOpt Print Level
%             0 = Silent
%             1 = Display one summary row
%             2 = Display one row for each unique global minima found. 
%             3 = Display one row for each unique global and local minima found.
%             4 = Display one row for each unique global and local minima and each solver failure.
%             5 = Display one row for each of the solver calls.
%             For PriLev 2,3,4 and 5 each row prints the following information:
%                 The minima are sorted, lowest value first (possibly the
%                 global minima)
%                 1.  #   ..    The trial number
%                 2.      ..    A code that classifies the solution 
%                               =1  Global minima, all with the same f(x) value (tolerance eps_f)
%                               =2,3,4,... The 2nd best, 3rd best, etc. group of local minima
%                               =-1,-2,... Feasible x_0, but solver failed to find feasible x_k
%                                         -1 has lowest f(x_0), -2 2nd best, etc.
%                               =-93 Solver failures, ranked from lowest f(x_k) + deviation
%                               =-92 Nonunique global and local minima
%                               =-91 Nonunique solver failure solutions
%                 3.  Ex  ..    ExitFlag returned from localSolver (normally 0)
%                 4.  Inf ..    Inform value returned from localSolver (normally 0)
%                 5.  f_x ..    Function value f(x_k) at final solution point x_k
%                 6.  f_0 ..    Function value f(x_0) at initial point x_0 minima
%                 7.  |x0-x| .. Distance from initial point x_0 to final point x_k, ||x_0-x_k||_2
%                 8.  NV  ..    Sum of deviation outside bounds in nonlinear constraints at x_k 
%                 9.  N0  ..    Sum of deviation outside bounds in nonlinear constraints at x_0
%                 10. LV  ..    Sum of deviation outside bounds in linear constraints at x_k
%                 11. L0  ..    Sum of deviation outside bounds in linear constraints at x_0
%                 8,9,10,11 are not displayed if NV or N0 <= cTol, or LV or L0 <= bTol
%             6 = PriLev 6 and higher rows with x_k is displayed
%             7 = PriLev 7 and higher rows with x_0 is displayed - before the x_k rows!
%             8 = tomRun (PrintResult) output from every optimization,
%                 print level 1
%             9 = tomRun (PrintResult) output from every optimization,
%                 print level 2
%                 For constrained problems output of sum(|constr|) and
%                 information if optimal point was accepted w.r.t.
%                 feasibility
%            10 = tomRun (PrintResult) output from every optimization,
%                 print level 3
%            11 = PriLevOpt = 1 are sent to the subsolver (GO.localSolver)
%            12 = PriLevOpt = 2 are sent to the subsolver
%                 In general max(0,PriLevOpt-10) are sent to the subsolver
%
% GO
%   localSolver The local solver used to run all local optimizations.
%               Default is the license dependent output of
%               GetSolver('con',1,0).
%
% optParam    Structure in Prob, Prob.optParam. Defines optimization
%             parameters. Fields used:
%   fGoal     Goal for function value f(x), if empty not used
%             If fGoal is reached, no further local optimizations are done
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal ~= 0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%             Default 1e-8.
%   bTol      Linear constraint feasibility tolerance. Default 1e-6
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:
%             If empty, all variables are assumed non-integer
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer
%             variables
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
% ExitFlag Exit flag 
%          = 0 normal output (Inform=0), of if fGoal set and found (Inform=8)
%          = 1 A solution reaching the user defined fGoal was not found (Inform=7)
%          = 4 Infeasible problem, or complete solver failure, see Inform 4,5,6
% Inform   Code telling type of convergence
%          = 0 Normal 
%          = 4 Linear infeasible problem
%          = 5 Nonlinear infeasible problem
%          = 6 Solver failure. Initial x_0 feasible
%          = 7 fGoal not reached 
%          = 8 fGoal reached
%          = 9 No initial x_0 possible to use. No optimization done
% DualGap  Relative duality gap, max(0,fIPMin-fLB) / fIPMin, if fIPMin ~=0
%            max(0,fIPMin-fLB) if fIPMin == 0. If fIPMin ~=0:
%            Scale with 100, 100*DualGap, to get the percentage duality gap.
%            Absolute value duality gap: scale with fIPMin, fIPMin * DualGap
% x_k      All global minima found
% f_k      Function value at global optimum
% c_k      Constraint values for all x_k (not computed)
% Solver   multiMin with local solver "localSolver"
% SolverAlgorithm  Description of method used
% Iter     Total number of iterations
% FuncEv   Total number of function evaluations
% GradEv   Total number of gradient evaluations
% HessEv   Total number of Hessian evaluations
% ConstrEv Total number of constraint function evaluations
%
% ---------------------------------------------------------------------
% The special field Result.multiMin is also used for warm start.
% Fields in Result.multiMin:
% ---------------------------------------------------------------------
%   The subfield Result.multiMin.Info contains counters etc:
%   Result.multiMin.Info:
%   localTry    = localTry+nTrial;
%   Iter        Total number of iterations
%   FuncEv      Total number of function evaluations
%   GradEv      Total number of gradient evaluations
%   HessEv      Total number of Hessian evaluations
%   ConstrEv    Total number of constraint function evaluations
%   ExitFlag    ExitFlag (converted to Result.Inform at return)
%   nGlobal     Number of global optima (same f(x) and feasible) found
%   nLocal      Number of local optima found (nGlobal + all other local minima)
%   nFail       Total number of solver failures
%   nSolution   Number of unique solutions with infeasible results
%   nFailF0OK   Number of infeasible solutions, where initial x_0 is feasible
%   nTrial      Total number of nonlinear problems solved, with different x_0
%   M           Total number of x_0 tried, M >= nTrial. Some x_0 rejected, if M > nTrial
%   localSolver The local solver used in all local optimization runs
%
%   Besides Result.multiMin.Info, the following fields are defined in Result.multiMin:
%   P           The trial numbers
%   LM          A code classifying the result of each nTrial optimization
%   X0          x_0 at each nTrial
%   XX          x_k found at each nTrial optimization
%   F0          f(x_0) at each nTrial. Row 2 f(x_0) + sum of linear and nonlinear infeasibilities
%   FX          f(x_k) at each nTrial. Row 2 f(x_k) + sum of linear and nonlinear infeasibilities
%   C0          Row 1 sum of linear, Row 2 sum of nonlinear infeasibility (at x_0) for each nTrial.
%   CX          Row 1 sum of linear, Row 2 sum of nonlinear infeasibility (at x_k) for each nTrial.
%   EX          Row 1 Code classifying solution, Row 2: ExitFlag, Row 3: Inform (from localSolve)
%               Code for each solver solution:
%               0 = not computed, f=NaN. 1 = non-unique infeasible solution, 
%               2 = non-unique feasible solution, 3 = infeasible solution, 
%               4 = infeasible solution, but feasible x_0, 5 = unique minima
%
% USAGE:
%
% Driver call, including printing with level 2:
%    Result = tomRun('multiMin',Prob,2);
%
% Direct solver call, Prob.xInit not set:
%    Result = multiMin(Prob);
%    generates min(3000,30*Prob.N) random initial points
% or with Prob.xInit set to 100, to generate 100 random initial points
%    Prob.xInit = 100;
%    Result = multiMin(Prob);
% or with a matrix of 202 generated initial points in [0,1]^Prob.N as input:
%    Prob.xInit = rand(Prob.N,202);
%    Result = multiMin(Prob);
% or with Prob.xInit set to -100, to generate 100 initial points
%    with Latin Hypercube random design
%
% Then print the result of the best local search with the call
%    PrintResult(Result);
%
% To make 20 repeated set of calls using warm start
%
%    Prob.xInit     = 100; % Use 100 initial points every call
%    Prob.PriLevOpt = 2;   % Display the list of minima last in multiMin
%    Result = tomRun('multiMin',Prob,2);
%    Prob.WarmStart = 1;
%    for i = 1:20
%        Prob   = WarmDefGLOBAL('multiMin',Prob,Result);
%        Result = tomRun('multiMin',Prob,2);
%        pause
%    end
%
% To run both a given set of initial values X, and 200 random initial values do
%    Prob.xInit     = X;   % Use initial points given in matrix X
%    Prob.PriLevOpt = 2;   % Display the list of minima last in multiMin
%    Result         = tomRun('multiMin',Prob,2);
%    Prob.WarmStart = 1;
%    % Note! Because of the warm start, the random generator seed will not be
%    % reset using Prob.RandState, it keeps the value from the calling routine
%    % If calling first with xInit = 200, then do warm start with xInit = X,
%    % then the random seed generator seed will be reset using Prob.RandState
%    Prob.xInit     = 200; % Use 200 initial points
%    Prob           = WarmDefGLOBAL('multiMin',Prob,Result);
%    Prob.Mx0       = 50;  % Select the 50 "best" of the 200 initial points
%    Result         = tomRun('multiMin',Prob,2);

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 13, 2006.    Last modified July 23, 2011.

function Result = multiMin(Prob, xInit)

if nargin < 2
   xInit = [];
   if nargin < 1
      error('multiMin needs input structure Prob');
   end
end

global NARG NLP_f NLP_c NLP_g NLP_dc NLP_r NLP_J NLP_x

n              = Prob.N;
probType       = Prob.probType;
CLS            = 0;
dLin           = Prob.mLin;
dNoL           = Prob.mNonLin;
if isfield(Prob.GO,'localSolver')
   localSolver = deblank(Prob.GO.localSolver);
else
   localSolver = [];
end
if isempty(localSolver) & probType==checkType('cls')
   SolvType    = checkType('cls');
   localSolver = GetSolver('cls',0,0);
   CLS         = 1;
else
   SolvType    = checkType('con');
end
if isempty(localSolver)
   % Always try to use snopt, not npsol
   if n > 200
      if checkLicense('snopt')
         localSolver = 'snopt';
      elseif checkLicense('npsol')
         localSolver = 'npsol';
      end
   else
      if checkLicense('snopt')
         localSolver = 'snopt';
      elseif checkLicense('npsol')
         localSolver = 'npsol';
      end
   end
   if isempty(localSolver)
      if checkLicense('knitro')
         localSolver = 'knitro';
      elseif checkLicense('conopt')
         localSolver = 'conopt';
      else
         localSolver = 'conSolve';
      end
   end
end
if checkLicense('minos')
   SolverQP    = 'qp-minos';
elseif checkLicense('knitro')
   SolverQP    = 'knitro';
elseif checkLicense('conopt')
   SolverQP    = 'conopt';
else
   SolverQP    = 'qld';
end
Prob = ProbCheck(Prob,localSolver,SolvType);
Prob = iniSolve(Prob,SolvType,1,1);

if strcmpi(localSolver,'snopt')
   SOLSNOPT    = 1;
else
   SOLSNOPT    = 0;
end
if strcmpi(localSolver,'npsol')
   SOLNPSOL    = 1;
else
   SOLNPSOL    = 0;
end

if isempty(Prob.FUNCS.g)
   DFREE       = 1;
else
   DFREE       = 0;
end

if DFREE
   MaxFunc     = DefPar(Prob.optParam,'MaxFunc',50000);
   xInit       = DefPar(Prob,'xInit',min(MaxFunc/25,30*n));
   if SOLSNOPT | SOLNPSOL
      Prob.SOL.optPar(9)  = 1E-4;
      Prob.SOL.optPar(10) = 1E-4;
      Prob.SOL.optPar(40) = 0;
      Prob.NumDiff        = 6;
      if dNoL > 0
         Prob.ConsDiff    = 6;
      end
      %Prob.SOL.optPar(22) = 0.1;
      %Prob.SOL.PrintFile  = 'snopt.txt';
      %Prob.SOL.optPar(1)  = 111111;
      %Prob.SOL.optPar(5)  = 1;
      %Prob.SOL.optPar(39) = 0;
   end
else
   MaxFunc     = Inf;
   xInit       = DefPar(Prob,'xInit',min(3000,30*n));
end
if any(size(xInit)) > 1
   % Check dimensions of given matrix
   if ~(size(xInit,1) == n)
      fprintf('Size of input matrix xInit is %d x %d\n',size(xInit));
      fprintf('Number of rows must be problem dimension %d\n',Prob.N);
      fprintf('Number of columns are equal to number of initial points\n');
      error('First dimension must have length == Prob.N');
   end
end
if length(xInit) <= 1
    SCALAR = 1;
else
    SCALAR = 0;
end

PriLev         = Prob.PriLevOpt;     % Print level
Prob.PriLevOpt = max(0,PriLev-10);   % Output from subsolver if PriLevOpt > 10
xEqTol         = DefPar(Prob,'xEqTol',1E-5);
fEqTol         = DefPar(Prob,'fEqTol',1E-7);
fCut           = DefPar(Prob,'fCut',1E10);
fGoal          = DefPar(Prob.optParam,'fGoal',-inf);
fTol           = Prob.optParam.eps_f;     % Relative tolerance for fGoal
x_L            = Prob.x_L(:);   % Lower bounds
x_U            = Prob.x_U(:);   % Upper bounds
x_LL           = x_L;
x_UU           = x_U;
c_L            = Prob.c_L(:);   % Lower bounds on nonlinear constraints
c_U            = Prob.c_U(:);   % Upper bounds on nonlinear constraints
b_L            = Prob.b_L(:);   % Lower bounds on linear constraints
b_U            = Prob.b_U(:);   % Upper bounds on linear constraints
K              = dLin + dNoL;
bTol           = Prob.optParam.bTol;      % Linear Constraint feasibility tolerance
if dLin > 0 & dNoL == 0
   cTol        = bTol;
else
   cTol        = Prob.optParam.cTol;
end

% Integer variables
IntVars         = DefPar(Prob.MIP,'IntVars',[]);
% Logical vector for integers
IV              = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('multiMin: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars         = find(IV);
nI              = length(IntVars);
RealVars        = find(~IV);
nR              = length(RealVars);
WarmStart       = DefPar(Prob,'WarmStart',0);

if nI == n
   nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
   fprintf('multiMin: Do not call multiMin for PURE integer programs\n');
   fprintf('Number of combinations = %d. ',nMax);
   if nMax < 100000
      fprintf('Call e.g. glcFast or glcDirect instead. ');
   else
      fprintf('Call e.g. minlpSolve instead. ');
   end
   fprintf('\n');
   error('multiMin: Wrong type of optimization problem');
end

if WarmStart == 1 & isfield(Prob,'multiMin')
   Info         = Prob.multiMin.Info;
   localTry     = DefPar(Info,'localTry',0);
   totIter      = DefPar(Info,'Iter',0);
   totfEval     = DefPar(Info,'FuncEv',0);
   gradEv       = DefPar(Info,'gradEv',0);
   HessEv       = DefPar(Info,'HessEv',0);
   totcEval     = DefPar(Info,'ConstrEv',0);
   nGlobal0     = DefPar(Info,'nGlobal',0);
   nLocal0      = DefPar(Info,'nLocal',0);
   nFail0       = DefPar(Info,'nFail',0);
   nSolution0   = DefPar(Info,'nSolution',0);
   nFailF0OK0   = DefPar(Info,'nFailF0OK',0);
   nTrial0      = DefPar(Info,'nTrial',0);
   M0           = DefPar(Info,'M',0);
   localSolver  = DefPar(Info,'localSolver',[]);

   if localTry == 0 | nTrial0 == 0
      WarmStart = 0;
   else
      %P0        = Prob.multiMin.P;
      LM0       = Prob.multiMin.LM;
      X00       = Prob.multiMin.X0;
      XX0       = Prob.multiMin.XX;
      %F00       = Prob.multiMin.F0;
      FX0       = Prob.multiMin.FX;
      %C00       = Prob.multiMin.C0;
      %CX0       = Prob.multiMin.CX;
      %EX0       = Prob.multiMin.EX;
   end
else
   Prob.WarmStart = 0;  % Must set 0 to avoid warm start in NLP solver
   WarmStart      = 0;
end
if WarmStart == 0
   localTry     = 0;
   totIter      = 0;
   totfEval     = 0;
   gradEv       = 0;
   HessEv       = 0;
   totcEval     = 0;

   if SCALAR
      if isfield(Prob,'RandState')
         RandState = Prob.RandState;
      else
         RandState = -1;
      end
      if isempty(RandState), RandState = -1; end
      % Set pseudo random generator if RandState ~=0
      if length(RandState) >1 | RandState > 0
         rand('state',RandState);
      elseif RandState <0
         rand('state',sum(100*clock));
      end
      % Print RandState and rand('state') 
      % fprintf('RandState %f sumState %f \n',sum(RandState),sum(rand('state')))
   end
end

Prob.WarmStart = 0;  % Must set 0 to avoid warm start in NLP solver

if SCALAR
    % Check for Inf and set to lower values.
    % Find lowest/highest finite bound
    xMin = min(x_L(isfinite(x_L)));
    xMax = max(x_U(isfinite(x_U)));
    if CLS
        L = 100;
    else
        L = 1000;
    end
    if isempty(xMin)
        x_L(isinf(x_L)) = min(-L,x_U(isinf(x_L))-L);
    else
        x_L(isinf(x_L)) = min(xMin,min(-L,x_U(isinf(x_L))-L));
    end
    if isempty(xMax)
        x_U(isinf(x_U)) = max(L,x_L(isinf(x_U))+L);
    else
        x_U(isinf(x_U)) = max(xMax,max(L,x_L(isinf(x_U))+L));
    end
    x_D = x_U-x_L;
    if xInit > 0
       M = xInit;
    else
       M = abs(xInit);
       xInit = daceInit(M,min(100,ceil(M/2)),x_L,x_U,[],0.05*norm(x_D));
       SCALAR = 0;
    end
    if nI > 0
       nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
       if nMax > M
          fprintf('multiMin: Do not call multiMin for MINLP with many IP combinations\n');
          fprintf('M = %d < ',M);
          fprintf('Number of combinations = %d. ',nMax);
          fprintf('Call e.g. multiMINLP instead. ');
          fprintf('\n');
          fprintf('            ');
          error('multiMin: Wrong type of optimization problem');
       end
    end
else
    M = size(xInit,2);
end

ExitFlag = 0;
absviol  = 0;
vTol     = 1E-8;
X0       = zeros(n,M);
XX       = zeros(n,M);
F0       = zeros(2,M);
FX       = zeros(2,M); 
C0       = zeros(2,M);
CX       = zeros(2,M);
EX       = zeros(3,M);
EX(2,:)  = NaN;
if WarmStart == 1
   Mx0   = DefPar(Prob,'Mx0',ceil(M/10));
   xD    = Inf*ones(M,1);      % Save distance to nearest old trial (x_0 or x_k)
end

if SCALAR
   if dLin > 0
      ProbLD = qpAssign(speye(n), zeros(n,1), Prob.A, b_L, b_U, ...
                        x_L, x_U, zeros(n,1), 'x0Feas');
      ProbLD.optParam = optParamDef('qp-minos',checkType('qp'));
      ProbLD = ProbCheck(ProbLD,'qp-minos',checkType('qp'));
   end
   if isempty(Prob.x_0)
      k1        = 0;
   else
      X0(:,1)   = Prob.x_0;
      k1        = 1;
      if dLin > 0
         [LV0 LVi] = consViol(Prob.A*Prob.x_0,b_L,b_U,absviol);
         C0(1,1)   = LV0;
      end
   end
   if ~isempty(IntVars)
      % nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
      if nMax <= M
         if nR == 0, M=nMax; end 
         blk = 0;
         % Generate all integer combinations
         IC = x_D(IntVars)+1 ;
         for j=1:nI
             i                  = IntVars(j);
             ic                 = IC(j);
             if j == 1
                blk             = ic;
                X0(i,k1+(1:ic)) = x_L(i) + (0:ic-1);
             else
                X0(i,k1+1:blk)     = x_L(i);
                for k = 2:ic
                    ix                             = k1+(k-1)*blk;
                    X0(IntVars(1:j-1),ix+1:ix+blk) = X0(IntVars(1:j-1),k1+[1:blk]);
                    X0(i             ,ix+1:ix+blk) = x_L(i) + k-1;
                end
                blk = blk*ic;
             end
         end
         % Find constraints in A only dependent on integer variables
         if dLin > 0
            ixI = find(all(Prob.A(:,RealVars)'==0));
            nixI = length(ixI);
            if nixI > 0
               % Check and delete any infeasible integer combinations
               LVi = consViolM(Prob.A(ixI,IntVars)*X0(IntVars,k1+(1:blk)),b_L(ixI),b_U(ixI),absviol);
               nDel = sum(LVi > 0);
               if nDel > 0
                  if PriLev > 1
                     fprintf('Deleted %d integer combinations out of %d\n',nDel,blk);
                  end
                  l    = blk-nDel;
                  % Move columns to be saved to the left
                  X0(IntVars,k1+(1:l)) = X0(IntVars,k1+find(LVi == 0));
                  blk  = l;
               end
            end
         end
         if blk==0
             error('Found no integer-feasible starting point.');
         end
         % Possibly shrinked blk.
         j = floor((M-k1)/blk);
         % IP block copied j times
         for k=2:j
             ix                      = k1+(k-1)*blk;
             X0(IntVars,ix+1:ix+blk) = X0(IntVars,k1+[1:blk]);
         end
         % Generate random IP points for the rest of the initial points
         % Simple solution, just take any of the blk combinations each time
         l = k1+blk*j;
         j = M-l;
         if j > 0
            ix                = max(1,ceil(blk*rand(j,1)));
            X0(IntVars,l+1:M) = X0(IntVars,k1+ix);
         end
      end
   end
   if nR > 0
      % Check if bounds huge, or squeezed
      ix = x_D(RealVars) >= 1000; 
      if any(ix) 
         HUGE = 1;
         y    = zeros(nR,1); % Define work area
      else
         HUGE = 0;
      end
      if HUGE
         % Find a base point close to 0
	 ixH           = find(ix);
	 ixI           = find(~ix);
	 xB            = zeros(nR,1);
	 if nI > 0
            xLL        = x_L(RealVars);
            xUU        = x_U(RealVars);
	 else
            xLL        = x_L;
            xUU        = x_U;
         end
         xB(ixH)       = max(xLL(ixH),min(xUU(ixH),xB(ixH)));
         xB(ixI)       = 0.5*(xLL(ixI)+xUU(ixI));
	 xSg           = zeros(nR,1);
	 xSg(xB==xLL)  = 1;
	 xSg(xB==xUU)  = -1;
	 xSg(ixI)      = 2;
	 ix            = find(xSg==0);
	 for j=1:length(ix)
             i         = ix(j);
	     if abs(xLL(i)) <= abs(xUU(i)/100)
	        xSg(i) = 1;
	        xB(i)  = xLL(i);
	     elseif abs(xUU(i)) <= abs(xLL(i))/100
	        xSg(i) = -1;
	        xB(i)  = xUU(i);
             end
         end
         xDD           = xUU-xLL;
	 ix            = find(xSg==0);
         xDD(ix)       = min(xUU(ix)-xB(ix),xB(ix)-xLL(ix));
	 Pow           = ceil(log10(max(xDD(ixH))+1));
	 Fac           = 10^(-Pow-1);
	 z             = zeros(nR,1);
	 ixHm          = find(xSg==-1);
	 ixHp          = find(xSg==1);
      end
   end
   for P = k1+1:M
       if nR > 0
          if HUGE > 0 
             y(ixI)          = (rand(length(ixI),1)-0.5).*xDD(ixI);
	     Fac             = 10*Fac;
	     if Fac > 1, Fac = 10^(-Pow); end
             y(ixH)          = Fac*randn(length(ixH),1).*xDD(ixH);
             y(ixHp)         = abs(y(ixHp));
             y(ixHm)         = -abs(y(ixHm));
             X0(RealVars,P)  = xB + y;
          else
             X0(RealVars,P) = x_L(RealVars) + rand(nR,1).*x_D(RealVars);
          end
       end
       if dLin > 0 
          z               = X0(:,P);
          [LV0 LVi]       = consViol(Prob.A*z,b_L,b_U,absviol);
          if LV0 > 0
             %fprintf('P %d LV0 %20.10f\n',P,LV0)
             ProbLD.QP.c  = -2*z;
             ProbLD.x_0   = z;
             ProbLD.P     = P;
             if ~isempty(IntVars)          % Fix integer variables
                ProbLD.x_L(IntVars) = z(IntVars);
                ProbLD.x_U(IntVars) = z(IntVars);
             end
             rLD          = tomRunFast(SolverQP,ProbLD);
             %PrintResult(rLD,2);
             z            = rLD.x_k;
             [LV0 LVi]    = consViol(Prob.A*z,b_L,b_U,absviol);
             C0(1,P)      = LV0;
             %if LV0 > 1E-10
             %   % PrintResult(rLD,2);
             %   fprintf('multiMin: QP Problem %d Linear infeasible, LV0 = %20.10f\n',P,LV0)
             %   % keyboard
             %end
          end
       end
   end
else
   for P=1:M
       z = xInit(:,P);
       if ~isempty(IntVars)                 % Assure integer variables are integer valued
          z(IntVars) = max(x_L(IntVars),min(x_U(IntVars),round(z(IntVars))));
       end
       if dLin > 0
         [LV0 LVi] = consViol(Prob.A*z,b_L,b_U,absviol);
          % Right now make now adjustment if linear infeasible
          C0(1,P)   = LV0;
       end
       X0(:,P) = z;
   end
   k1 = 0;
end
remove = 0;
if WarmStart == 1
   for P = 1:M
       D = checkDist(X0(:,P),X00);
       if D < xEqTol
          remove  = remove + 1;
          F0(1,P) = NaN;
       else
          D = min(D,checkDist(X0(:,P),XX0));
          if D < xEqTol
             remove  = remove + 1;
             F0(1,P) = NaN;
          end
      end
      xD(P) = D;
   end
end

NARG  = []; NLP_x = []; NLP_f = []; NLP_c = []; NLP_g = []; NLP_dc = []; NLP_r = []; NLP_J = []; 
if Prob.Threads > 1
   F0(1,:) = multiMin_1_par(X0,Prob);
else
   for P = 1:M
       F0(1,P) = nlp_f(X0(:,P), Prob);         % Compute f(x_0) for each x_0
   end
end

if M > 1
   F0(1,k1+find(F0(1,k1+1:M) >= fCut)) = NaN;                    % Remove all f0 > fCut
end

% Remove same/similar initial values
remove = 0;
for P = 1:M-1
    f0   = F0(1,P);
    if ~isnan(f0)
       ix   = P+find(abs(f0-F0(1,P+1:M)) < fEqTol*max(1,abs(f0)));
       if ~isempty(ix)
          iE              = checkEq(X0(:,P),X0(:,ix),xEqTol);
          if ~isempty(iE)
             remove       = remove + length(iE);
             F0(1,ix(iE)) = NaN;
          end
       end
    end
end

if WarmStart == 1
   ix              = find(~isnan(F0(1,:)));
   [Dmax iMax]     = sort(xD(ix),'descend');
   if Mx0 < length(ix)
      F0(1,ix(iMax(Mx0+1:end))) = NaN;        % Remove x_0 closest to old x_0 and x_k
   end
   if PriLev > 1
      fprintf('Keep Mx0=%d. ',Mx0);
      fprintf('Delete %d x_0 closest to old initial and final points\n',length(ix)-Mx0);
   end
end
if k1 == 0 | (k1 == 1 & F0(1,1) >= fCut)
   ix                 = find(~isnan(F0(1,:)));
   if isempty(ix)
      ix              = 1;
   else
      [F00,ixS]       = sort(F0(1,ix));
      ix              = ix(ixS);
   end
else
   ix                 = find(~isnan(F0(1,2:end)));
   if isempty(ix)
      ix              = 1;
   else
      ix              = ix+1;
      [F00,ixS]       = sort(F0(1,ix));
      ix              = [1,ix(ixS)];    % Given x_0 used in trial 1
      if F0(1,1) > F0(1,ix(2))
         ix(1)           = ix(2);       % Given x_0 used in trial 2
         ix(2)           = 1;
      end
   end
end
Name               = Prob.Name;
NV                 = 0;
LV                 = 0;

if Prob.Threads > 1
   [Eval,C02,CX1,CX2,FX1,EX2,EX3,XX1] = multiMin_2_par(Prob,X0(:,ix),F0(:,ix), ...
         IntVars,localSolver,SOLSNOPT,PriLev);
   totfEval   = totfEval + Eval(1);
   gradEv     = gradEv   + Eval(2);
   totcEval   = totcEval + Eval(3);
   HessEv     = HessEv   + Eval(4);
   totIter    = totIter  + Eval(5);
   C0(2,ix)   = C02;
   CX(1,ix)   = CX1;
   CX(2,ix)   = CX2;
   EX(2,ix)   = EX2;
   EX(3,ix)   = EX3;
   FX(1,ix)   = FX1;
   XX(:,ix)   = XX1;
else
   % Local Optimization for each x_0
   nix = length(ix);
   for i=1:nix
       P     = ix(i);
       x_0   = X0(:,P);
       NLP_x = x_0;
       NLP_f = F0(1,P);
       if dNoL > 0
          c_0       = nlp_c(x_0, Prob);
          [NV0 NVi] = consViol(c_0,c_L,c_U,absviol);
          %C0(2,P)   = NV0;
          C02(i)   = NV0;
       end
       Prob.x_0     = x_0;
       Prob.Name    = [Name, ' - Trial ' num2str(P)];
       if ~isempty(IntVars)          % Fix integer variables
          Prob.x_L(IntVars) = x_0(IntVars);
          Prob.x_U(IntVars) = x_0(IntVars);
       end
       r            = tomRun(localSolver,Prob,PriLev-7);
       %PrintResult(r,2);
       if SOLSNOPT == 1 & r.FuncEv == 0  % NOT NEEDED: & Inform ~= -999
          totfEval  = totfEval + r.iwCount(3);
          gradEv    = gradEv   + sum(r.iwCount(4:6));
          totcEval  = totcEval + r.iwCount(7);
       else
          totfEval  = totfEval + r.FuncEv;
          gradEv    = gradEv   + r.GradEv;
          totcEval  = totcEval + r.ConstrEv;
       end
       HessEv       = HessEv   + r.HessEv;
       totIter      = totIter  + r.Iter;
       Inform       = r.Inform;
       f_k          = r.f_k;
       FX(1,P)      = f_k;
       EX(2,P)      = r.ExitFlag;
       EX(3,P)      = Inform;
       % Avoid tiny value outside bounds
       x_k          = min(Prob.x_U,max(Prob.x_L,r.x_k));
       XX(:,P)      = x_k;
       if dLin > 0                   % Linear constraint violation
          z         = r.Ax;
          if isempty(z)
             z      = Prob.A*x_k;
          end
          [LV LVi]  = consViol(z,b_L,b_U,absviol);
          CX(1,P)   = LV;
       end
       if dNoL > 0                  % Nonlinear constraint violation
          c_k       = r.c_k;
          if isempty(c_k)
             c_k    = nlp_c(x_k, Prob);
          end
          [NV NVi]  = consViol(c_k,c_L,c_U,absviol);
          CX(2,P)   = NV;
       end
       if ~isinf(fGoal) & LV <= dLin*bTol & NV <= dNoL*cTol % Check if goal is fulfilled
          if fGoal == 0
             if abs(f_k-fGoal) <= fTol
                 ExitFlag = 1;
                 break;
             end
          elseif abs(f_k-fGoal) <= abs(fGoal) * fTol
               ExitFlag = 1;
             break;
          end
       end
       if (SOLSNOPT & any(Inform == [11 12])) | (SOLNPSOL & Inform == 2)
          if P > 1  % Check that previous problems not solved OK
             FAIL =  all(sum(CX(:,1:P-1),1) > vTol);
          else
             FAIL =  1;
          end
          if all(C0(1,:) >= vTol) & FAIL
             % Linear constraints infeasible - waste of time to continue
             F0(1,P+1:end) = NaN;   % Mark the rest of the problems not solved
             ExitFlag      = 4;
             break
          end
       end
       if totfEval > MaxFunc
          F0(1,ix(i+1:end)) = NaN;
          break; 
       end
   end
end
% Reset lower/upper bounds for integer variables
Prob.x_L(IntVars) = x_L(IntVars);
Prob.x_U(IntVars) = x_U(IntVars);

% Classify all minima: 
%    0 = not computed, f=NaN. 1 = non-unique infeasible solution,        2 = non-unique feasible solution
%    3 = infeasible solution, 4 = infeasible solution, but feasible x_0, 5 = unique minima

vTol                         = K*cTol;
z                            = sum(CX,1);
v                            = sum(C0,1);
F0(2,:)                      = F0(1,:)+v;
FX(2,:)                      = FX(1,:)+z;
EX(1,z <= vTol)              = 5; 
EX(1,z >  vTol)              = 3; 
EX(1,isnan(F0(1,:)))         = 0; 
EX(1,z >  vTol & v <= vTol)  = 4; 


if WarmStart == 1
   % Remove any minima that is equal to a minimum from last run
   remove = 0;
   ixO    = find(EX(1,:)==5);
   for j=1:nLocal0
       fX = FX0(2,j);
       ix = find(abs(fX-FX(2,ixO)) < fEqTol*max(1,abs(fX)));
       if ~isempty(ix)
          iE   = checkEq(XX0(:,j),XX(:,ixO(ix)),xEqTol);
          if ~isempty(iE)
             remove              = remove + length(iE);
             P       = ixO(ix(iE));
             EX(1,P) = 2;
             ixO     = setdiff(ixO,P);
          end
       end
   end

   % Remove any failed solution that is equal to a failure from last run
   remove = 0;
   ixO    = find(EX(1,:)==3);
   for j=nLocal0+1:nSolution0
       fX = FX0(2,j);
       ix = find(abs(fX-FX(2,ixO)) < fEqTol*max(1,abs(fX)));
       if ~isempty(ix)
          iE   = checkEq(XX0(:,j),XX(:,ixO(ix)),xEqTol);
          if ~isempty(iE)
             remove              = remove + length(iE);
             P       = ixO(ix(iE));
             EX(1,P) = 1;
             ixO     = setdiff(ixO,P);
          end
       end
   end
end

% Find fOpt - sort minima, remove duplicates
ix             = find(EX(1,:)==5);
if ~isempty(ix)
   [fOpt2 ixO]  = sort(FX(2,ix));
   ixOpt        = ix(ixO);
   fOpt         = FX(1,ixOpt);
   remove       = 0;
   k            = length(fOpt);
   i            = 1;
   while i < k
      P        = ixOpt(i);
      fX       = fOpt(i);
      ix       = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
      if ~isempty(ix)
         iE    = checkEq(XX(:,P),XX(:,ixOpt(ix)),xEqTol);
         if ~isempty(iE)
            remove              = remove + length(iE);
            EX(1,ixOpt(ix(iE))) = 2;
         end
         i     = i + length(ix);
      end
      i        = i + 1;
   end
end
   
% Find infeasible solutions - remove duplicates
ix             = find(EX(1,:)==3);
if ~isempty(ix)
   [fOpt ixO]  = sort(FX(2,ix));
   ixOpt       = ix(ixO);
   remove      = 0;
   k           = length(fOpt);
   i           = 1;
   while i < k
       P       = ixOpt(i);
       fX      = FX(2,P);
       ix      = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
       if ~isempty(ix)
          iE   = checkEq(XX(:,P),XX(:,ixOpt(ix)),xEqTol);
          if ~isempty(iE)
             remove              = remove + length(iE);
             EX(1,ixOpt(ix(iE))) = 1;
          end
          i    = i + length(ix);
       end
       i       = i + 1;
   end
end

% Find fOpt - sort minima, classify
ix             = find(EX(1,:)==5);
idx            = [];
LM             = zeros(1,M);
if ~isempty(ix)
   [fOpt ixO]  = sort(FX(2,ix));
   ixOpt       = ix(ixO);
   k           = length(fOpt);
   if k == 1
      P        = ixOpt(1);
      idx      = [idx,P];
      LM(P)    = 1;
   else
      LMi          = 0;            % Global Opt = 1, other local minima 2,3,...
      i            = 1;
      while i < k
         LMi   = LMi+1;
         P     = ixOpt(i);
         idx   = [idx,P];
         LM(P) = LMi;
         fX    = FX(2,P);
         ix    = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
         if ~isempty(ix)
            idx           = [idx,ixOpt(ix)];
            LM(ixOpt(ix)) = LMi;
            i             = i + length(ix);
         end
         i     = i + 1;
         if i == k
            P      = ixOpt(k);
            idx    = [idx,P];
            LMi    = LMi+1;
            LM(P)  = LMi;
         end
      end
   end
end

% Find software failures, with feasible initial points and sort
ix             = find(EX(1,:)==4);
if ~isempty(ix)
   [fOpt ixO]  = sort(F0(1,ix));
   ixOpt       = ix(ixO);
   k           = length(fOpt);
   if k == 1
      P        = ixOpt(1);
      idx      = [idx,P];
      LM(P)    = -1;
   else
      LMi      = 0;
      i        = 1;
      while i < k
         LMi   = LMi-1;
         P     = ixOpt(i);
         idx   = [idx,P];
         LM(P) = LMi;
         fX    = F0(1,P);
         ix    = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
         if ~isempty(ix)
            idx           = [idx,ixOpt(ix)];
            LM(ixOpt(ix)) = LMi;
            i             = i + length(ix);
         end
         i     = i + 1;
         if i == k
            P      = ixOpt(k);
            idx    = [idx,P];
            LMi    = LMi-1;
            LM(P)  = LMi;
         end
      end
   end
end

% Find other solutions and sort
ix             = find(EX(1,:)==3);
if ~isempty(ix)
   [fOpt ixO]  = sort(FX(2,ix));
   ixOpt       = ix(ixO);
   k           = length(fOpt);
   if k == 1
      P        = ixOpt(1);
      idx      = [idx,P];
      LM(P)    = -93;
   else
      i            = 1;
      while i < k
         P     = ixOpt(i);
         idx   = [idx,P];
         LM(P) = -93;
         fX    = FX(2,P);
         ix    = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
         if ~isempty(ix)
            idx           = [idx,ixOpt(ix)];
            LM(ixOpt(ix)) = -93;
            i             = i + length(ix);
         end
         i     = i + 1;
         if i == k
            P      = ixOpt(k);
            idx    = [idx,P];
            LM(P)  = -93;
         end
      end
   end
end

% Find duplicate local solutions, and sort 
ix              = find(EX(1,:)==2);
if ~isempty(ix)
   [fOpt ixO]   = sort(FX(2,ix));
   ixOpt        = ix(ixO);
   k            = length(fOpt);
   if k == 1
      P         = ixOpt(1);
      idx       = [idx,P];
      LM(P)     = -92;
   else
      i         = 1;
      while i < k
         P      = ixOpt(i);
         idx    = [idx,P];
         LM(P)  = -92;
         fX     = FX(2,P);
         ix     = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
         if ~isempty(ix)
            idx           = [idx,ixOpt(ix)];
            LM(ixOpt(ix)) = -92;
            i             = i + length(ix);
         end
         i      = i + 1;
         if i == k
            P      = ixOpt(k);
            idx    = [idx,P];
            LM(P)  = -92;
         end
      end
   end
end

% Find duplicate nonfeasible solutions 
ix              = find(EX(1,:)==1);
if ~isempty(ix)
   [fOpt ixO]   = sort(FX(2,ix));
   ixOpt        = ix(ixO);
   k            = length(fOpt);
   if k == 1
      P         = ixOpt(1);
      idx       = [idx,P];
      LM(P)     = -91;
   else
      i         = 1;
      while i < k
         P      = ixOpt(i);
         idx    = [idx,P];
         LM(P)  = -91;
         fX     = FX(2,P);
         ix     = i+find(abs(fX-fOpt(i+1:k)) < fEqTol*max(1,abs(fX)));
         if ~isempty(ix)
            idx           = [idx,ixOpt(ix)];
            LM(ixOpt(ix)) = -91;
            i             = i + length(ix);
         end
         i      = i + 1;
         if i == k
            P      = ixOpt(k);
            idx    = [idx,P];
            LM(P)  = -91;
         end
      end
   end
end
if WarmStart == 1
   nLocal    = sum(LM > 0);
   fOpt      = [FX0(2,1:nLocal0),FX(2,idx(1:nLocal))];
   nTot      = length(fOpt);
   [fO iOpt] = sort(fOpt);
   Midx0     = find(iOpt <=nLocal0);
   Midx      = find(iOpt > nLocal0);
   % Classify the merged minima
   LMM       = zeros(1,nLocal+nLocal0);
   if size(LMM,2) > 0
      k           = length(fO);
      if k == 1
         LMM(1)        = 1;
      else
         LMi           = 0;            % Global Opt = 1, other local minima 2,3,...
         i             = 1;
         while i < k
            LMi        = LMi+1;
            LMM(i)     = LMi;
            fX         = fO(1);
            ix         = i+find(abs(fX-fO(i+1:k)) < fEqTol*max(1,abs(fX)));
            if ~isempty(ix)
               LMM(ix) = LMi;
               i       = i + length(ix);
            end
            i          = i + 1;
            if i == k
               LMi     = LMi+1;
               LMM(k)  = LMi;
            end
         end
      end
   end
   LM0(1:nLocal0)      = LMM(Midx0);
   LM(idx(1:nLocal))   = LMM(Midx);
   ix                  = find(LM0 < 0 & LM0 > -91);
   if ~isempty(ix)
      Midx0            = [Midx0,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM < 0 & LM > -91);
   if ~isempty(ix)
      Midx             = [Midx,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM0 == -93);
   if ~isempty(ix)
      Midx0            = [Midx0,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM == -93);
   if ~isempty(ix)
      Midx             = [Midx,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM0 == -92);
   if ~isempty(ix)
      Midx0            = [Midx0,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM == -92);
   if ~isempty(ix)
      Midx             = [Midx,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM0 == -91);
   if ~isempty(ix)
      Midx0            = [Midx0,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
   ix                  = find(LM == -91);
   if ~isempty(ix)
      Midx             = [Midx,nTot + [1:length(ix)]];
      nTot             = nTot + length(ix);
   end
end

if PriLev > 1
       %fprintf('Try%3d ',i)
   for P = idx
       if PriLev == 2 & LM(P) ~= 1
          break;
       end
       if PriLev == 3 & LM(P) < 0 
          break;
       end
       if PriLev == 4 & LM(P) < -90 
          break;
       end
       fprintf(' ')
       fprintf('#%3d. ',P)
       fprintf('%3d ',LM(P))
       fprintf('Ex%2d. ',EX(2,P))
       fprintf('Inf%3d. ',EX(3,P))
       fprintf('f_x %17.10f ',FX(1,P))
       fprintf('f_0 %17.10f. ',F0(1,P))
       fprintf('|x0-x|%12.6f. ',norm(X0(:,P)-XX(:,P)));
       if dNoL > 0
          v = CX(2,P);
          if v > cTol
             fprintf('NV%10.2e ',v);
          else
             fprintf('NV%10.2e ',v);
	     %HKH
             %fprintf('             ');
          end
          v = C0(2,P);
          if v > cTol
             fprintf('N0%10.2e ',v);
          else
             fprintf('             ');
          end
       end
       if dLin > 1E-14
          v = CX(1,P);
          if v > bTol
             fprintf('LV%10.2e ',v);
          else
             fprintf('LV%10.2e ',v);
	     %HKH
             %fprintf('             ');
          end
          v = C0(1,P);
          if v > bTol
             fprintf('L0%10.2e ',v);
          else
             fprintf('             ');
          end
       end
       fprintf('\n')
       if PriLev > 6
          xprint(X0(:,P),'x_0');
       end
       if PriLev > 5
          xprint(XX(:,P),'x_k');
       end
    end
end

Result                 = ResultDef(Prob);
Result.Solver          = ['multiMin with local solver ' localSolver];
Result.SolverAlgorithm = 'Find local optima using multistart local search';

nGlobal            = sum(LM == 1);
nLocal             = sum(LM > 0);
nFail              = sum(LM == -93) + sum(LM < 0 & LM >= -91);
nFailF0OK          = sum(LM == -1);
nSolution          = length(idx);
nTrial             = sum(~isnan(F0(1,:)));

if PriLev > 0
   fprintf('   multiMin: ');
   if WarmStart == 0
      fprintf('%d global minima, ',nGlobal);
      fprintf('%d minima. ',nLocal);
      fprintf('%d solved ',nSolution);
      fprintf('out of %d tries, ',nTrial);
      fprintf('%d failed. ',nFail);
      if nTrial < M
         fprintf('%d deleted. ',M-nTrial);
      end
      fprintf('TotFuncEv %d. ',totfEval)
      if totcEval > 0
         fprintf('TotcEval %d. ',totcEval)
      end
      fprintf('xEqTol %8.2e',xEqTol)
      fprintf('\n');
   else
      fprintf('%d(%d) global minima, ',nGlobal,nGlobal0);
      fprintf('%d(%d) minima. ',nLocal,nLocal0);
      fprintf('%d(%d) solved ',nSolution,nSolution0);
      fprintf('out of %d(%d) tries, ',nTrial,nTrial0);
      fprintf('%d(%d) failed. ',nFail,nFail0);
      fprintf('TotFuncEv %d. ',totfEval)
      if totcEval > 0
         fprintf('TotcEval %d. ',totcEval)
      end
      fprintf('xEqTol %8.2e',xEqTol)
      fprintf('\n');
   end
end
Result.FuncEv       = totfEval;
Result.GradEv       = gradEv;
Result.HessEv       = HessEv;
Result.Iter         = totIter;
Result.ConstrEv     = totcEval;

Info.localTry       = localTry+nTrial;
Info.Iter           = totIter;
Info.FuncEv         = totfEval;
Info.GradEv         = gradEv;
Info.HessEv         = HessEv;
Info.ConstrEv       = totcEval;

if WarmStart == 1
   Info.nGlobal     = nGlobal   + nGlobal0;
   Info.nLocal      = nLocal    + nLocal0;
   Info.nFail       = nFail     + nFail0;
   Info.nSolution   = nSolution + nSolution0;
   Info.nFailF0OK   = nFailF0OK + nFailF0OK0;
   Info.nTrial      = nTrial    + nTrial0;
   Info.M           = M         + M0;
else
   Info.nGlobal     = nGlobal;
   Info.nLocal      = nLocal;
   Info.nFail       = nFail;
   Info.nSolution   = nSolution;
   Info.nFailF0OK   = nFailF0OK;
   Info.nTrial      = nTrial;
   Info.M           = M;
end

Info.localSolver     = localSolver;
Result.multiMin.Info = Info;

if WarmStart == 1
   nTot                          = nTrial0 + nTrial;
   Result.multiMin.P             = zeros(1,nTot);
   Result.multiMin.LM            = zeros(1,nTot);
   Result.multiMin.X0            = zeros(n,nTot);
   Result.multiMin.XX            = zeros(n,nTot);
   Result.multiMin.F0            = zeros(2,nTot);
   Result.multiMin.FX            = zeros(2,nTot); 
   Result.multiMin.C0            = zeros(2,nTot);
   Result.multiMin.CX            = zeros(2,nTot);
   Result.multiMin.EX            = zeros(3,nTot);

   Result.multiMin.P(1,Midx)     = idx;
   Result.multiMin.LM(:,Midx)    = LM(:,idx);
   Result.multiMin.X0(:,Midx)    = X0(:,idx);
   Result.multiMin.XX(:,Midx)    = XX(:,idx);
   Result.multiMin.F0(:,Midx)    = F0(:,idx);
   Result.multiMin.FX(:,Midx)    = FX(:,idx);
   Result.multiMin.C0(:,Midx)    = C0(:,idx);
   Result.multiMin.CX(:,Midx)    = CX(:,idx);
   Result.multiMin.EX(:,Midx)    = EX(:,idx);

   Result.multiMin.P(1,Midx0)    = Prob.multiMin.P;
   Result.multiMin.LM(:,Midx0)   = LM0;
   Result.multiMin.X0(:,Midx0)   = X00;
   Result.multiMin.XX(:,Midx0)   = XX0;
   Result.multiMin.F0(:,Midx0)   = Prob.multiMin.F0;
   Result.multiMin.FX(:,Midx0)   = FX0;
   Result.multiMin.C0(:,Midx0)   = Prob.multiMin.C0;
   Result.multiMin.CX(:,Midx0)   = Prob.multiMin.CX;
   Result.multiMin.EX(:,Midx0)   = Prob.multiMin.EX;

   if Info.nGlobal > 0
      Result.x_k      = Result.multiMin.XX(:,1:Info.nGlobal);
      Result.f_k      = Result.multiMin.FX(1,1);
      ExitFlag        = 0;
   elseif Info.nFailF0OK > 0
      % Solver failed to find feasible solution
      % There are some feasible initial values, pick the best
      ix              = Info.nLocal+[1+Info.nFailF0OK];
      Result.x_k      = Result.multiMin.X0(:,ix);
      Result.f_k      = Result.multiMin.F0(1,ix(1));
      ExitFlag        = 6;
   elseif ExitFlag == 4
      Result.x_k      = [];
      Result.f_k      = NaN;
   else
      Result.x_k      = Result.multiMin.XX(:,1);
      Result.f_k      = Result.multiMin.FX(1,1);
      ExitFlag        = 5;
   end
else
   Result.multiMin.P             = idx;
   Result.multiMin.LM            = LM(:,idx);
   Result.multiMin.X0            = X0(:,idx);
   Result.multiMin.XX            = XX(:,idx);
   Result.multiMin.F0            = F0(:,idx);
   Result.multiMin.FX            = FX(:,idx);
   Result.multiMin.C0            = C0(:,idx);
   Result.multiMin.CX            = CX(:,idx);
   Result.multiMin.EX            = EX(:,idx);

   if nGlobal > 0
      Result.x_k      = XX(:,idx(1:nGlobal));
      Result.f_k      = FX(1,idx(1));
      ExitFlag        = 0;
   elseif nFailF0OK > 0
      % Solver failed to find feasible solution
      % There are some feasible initial values, pick the best
      Result.x_k      = X0(:,idx(1:nFailF0OK));
      Result.f_k      = F0(1,idx(1));
      ExitFlag        = 6;
   elseif ExitFlag == 4
      Result.x_k      = [];
      Result.f_k      = NaN;
   elseif nTrial == 0
      Result.x_k      = [];
      Result.f_k      = NaN;
      ExitFlag        = 999;
   else
      Result.x_k      = XX(:,idx(1));
      Result.f_k      = FX(1,idx(1));
      ExitFlag        = 5;
   end
end
x_k = Result.x_k;
nxk = size(x_k,2);
if dLin > 0 && ~isempty(x_k)
   Result.Ax        = Prob.A*x_k;
else
   Result.Ax = [];
end
if dNoL > 0
   Result.c_k       = zeros(dNoL,nxk);
   for i=1:nxk
       Result.c_k(:,i) = nlp_c(x_k(:,i),Prob);
   end
end

if PriLev > 13
   xprint(x_k,'x_k:')
   if dLin > 0
      disp('[b_L Result.Ax b_U]')
      disp([b_L Result.Ax b_U])
   end
   if dNoL > 0
      disp('[c_L Result.c_k c_U]')
      disp([c_L Result.c_k c_U])
   end
end


if ExitFlag == 1
    Result.ExitText = ...
        ['fGoal ' num2str(fGoal) ' found, ' num2str(nTrial) ' tries. ' ...
        'Found ' num2str(nGlobal), ' global, ', num2str(nLocal) ' minima'];
    Result.ExitFlag = 0;
    Result.Inform   = 8;
elseif ExitFlag == 4
    Result.ExitText = 'Linear infeasible problem';
    Result.ExitFlag = 4;
    Result.Inform   = 4;
elseif ExitFlag == 5
    Result.ExitText = 'Nonlinear infeasible problem';
    Result.ExitFlag = 4;
    Result.Inform   = 5;
elseif ExitFlag == 6
    Result.ExitText = 'Solver failure. Initial x_0 feasible';
    Result.ExitFlag = 4;
    Result.Inform   = 6;
elseif ExitFlag == 999
    Result.ExitText = 'No initial x_0 possible to use';
    Result.ExitFlag = 4;
    Result.Inform   = 9;
elseif ~isempty(fGoal) & ~isinf(fGoal)
    Result.ExitText = ...
        ['fGoal ' num2str(fGoal) ' NOT found, ' num2str(nTrial) ' tries. ' ...
        'Found ' num2str(nGlobal), ' global, ', num2str(nLocal) ' minima'];
    Result.ExitFlag = 1;
    Result.Inform   = 7;
elseif localTry > 0
    Result.ExitText = ...
        ['Did ' num2str(nTrial) ' new local tries, ' num2str(localTry) ' previous. ' ...
        'Found ' num2str(nGlobal) ' (' num2str(nGlobal0) ') global, ' ...
         num2str(nLocal) '(' num2str(nLocal0) ') minima'];
    Result.ExitFlag = ExitFlag;
    Result.Inform   = 0;
else
    Result.ExitText = ...
        ['Did ' num2str(nTrial) ' local tries. ' ...
        'Found ' num2str(nGlobal), ' global, ', num2str(nLocal) ' minima'];
    Result.ExitFlag = ExitFlag;
    Result.Inform   = 0;
end
Result.ExitText = [ Result.ExitText '. TotFuncEv ' num2str(totfEval)];
if Prob.mNonLin > 0
   Result.ExitText = [ Result.ExitText '. TotConstrEv ' num2str(totcEval)];
end

if localTry > 0
    multiMin.ExitText = ['Tried ' num2str(nTrial) ' local searches, ' num2str(localTry) ' previous. '...
        'Found ' num2str(nGlobal), ' global, ', num2str(nLocal) ' minima'];
else
    multiMin.ExitText = ['Tried ' num2str(nTrial) ' local searches. ' ...
        'Found ' num2str(nGlobal), ' global, ', num2str(nLocal) ' minima'];
end

%HKH fix
%fprintf('   multiMin - WarmStart %d ExitText: %s\n',WarmStart,Result.ExitText)

Result     = endSolve(Prob,Result);
Result.f_k = Result.f_k-Prob.fConstant;

% -----------------------------------------------
function ix = checkEq(x_k,X,xTol)
% -----------------------------------------------
% Borde nyttja eps_x, kolla varje komponent
D  = tomsol(30,x_k,X)/max(1,norm(x_k));
%xprint(D)
L = D <= xTol;
% Check if more solutions are possible to delete
for i=1:length(L)
    if L(i) == 0
       ix = find(L(i+1:end)==0 & (D(i)-D(i+1:end)) <= 1E-5*D(i));
       xN = max(1,norm(X(:,i)));
       for j=1:length(ix)
           k = i + ix(j);
           if tomsol(30,X(:,i),X(:,k)) <= xTol*xN
              L(k) = 1;
           end
       end
    end
end
ix = find(L);

% -----------------------------------------------
function D = checkDist(x_k,X)
% -----------------------------------------------
D  = min(tomsol(30,x_k,X)/max(1,norm(x_k)));

% -----------------------------------------------
function [h_L1, h] = consViol(z,L,U,absviol)
% -----------------------------------------------
if absviol == 1
   h =-min(0,z-L)-min(0,U-z);
else
   h =-min(0,(z-L)./max(1,abs(L)))-min(0,(U-z)./max(1,abs(U)));
end
h_L1 = sum(h);

%fprintf('h_L1  %25.10f ',h_L1);
%h_L2 = norm(h);
%fprintf('h_L2  %25.10f\n',h_L2);

% -----------------------------------------------
function h_L1 = consViolM(Z,L,U,absviol)
% -----------------------------------------------
m             = size(Z,2);
h_L1          = zeros(m,1);
if absviol == 1
   for i=1:m
      h_L1(i) = sum(-min(0,Z(:,i)-L)-min(0,U-Z(:,i)));
   end
else
   for i=1:m
      h_L1(i) = sum(-min(0,(Z(:,i)-L)./max(1,abs(L)))-min(0,(U-Z(:,i))./max(1,abs(U))));
   end
end

% MODIFICATION LOG:
%
% 060713  hkh  Written
% 060714  hkh  All information in one Result output, xInit also as Prob.xInit
% 060715  hkh  Test on 1st dimension of xInit, must be == Prob.N
% 060715  hkh  Add ExitFlag, ExitText, Inform information, and warm start
% 060716  hkh  Try 100 times to get linear feasible when random initial x_0
% 060716  hkh  Add RandState to enable setting of rand('state',.)
% 060716  hkh  Prob.x_0 was not set OK for integer variables
% 060717  hkh  Prob.xEqTol is used to check if x_k equal to xOpt
% 060818  hkh  Always set a default on fCut(=1E10) to avoid stops at huge f
% 060818  hkh  multiMin detects NLLS, and uses default CLS solver.
% 060818  hkh  Safeguarded settings of infinite x_L, x_U
% 060818  hkh  Always use given Prob.x_0 as 1st trial, instead of random point
% 060819  hkh  Changed isscalar to simple length test, seems it was 7.x feature
% 061123  hkh  Insert one extra PriLev, 1, for one single summary row of output
% 070222  hkh  Revise code and format for IntVars
% 070815  med  Moved SCALAR definition to make sure warm start works
% 070815  med  Correct help text
% 071003  hkh  Added comments about warm start
% 071003  hkh  Avoid tiny value outside bounds: x_k = min(x_U,max(x_L,r.x_k));
% 071005  hkh  Change default tries to min(3000,30*dim(x))
% 080513  med  Now possible to use multiMin with MAD
% 080518  med  Always return a proper result struct Result = r (if no solution)
% 080529  med  Corrected fConstant
% 080827  hkh  Safe guard cTol in case of only linear constraints, cTol undef
% 081104  hkh  Use ResultDef to define all fields in r, avoid crash in endSolve
% 090827  hkh  Add ExitFlag 4, infeasible problem
% 090827  hkh  Detect infeasibility in 1st local run (SNOPT). Quit directly 
% 090912  hkh  Compute nFuncEv,nConstrEv, nGradEv from iwCount (only SNOPT)
% 090912  hkh  Add TotFuncEv and TotConstrEv to ExitText
% 090914  hkh  |xInit| x_0 pnts with Latin Hypercube, if xInit < 0
% 090917  hkh  Limit # of intervals, s=min(100,M/2), in call to daceInit
% 090918  hkh  Set r.FuncEv=1 if f_0 > fCut is computed
% 091008  hkh  Stop if nonlinear infeasible for n+1 tries, with no feasible f(x)
% 091013  hkh  Prob.fEqTol is used to check if f_k values are the same
% 091014  hkh  Complete revision. Classify and sort all minima. New output/warm start handling
% 091014  hkh  Identify and remove equal minima. Find and return all global optima.
% 091018  hkh  Fix integer variables in QP solve. Reset lower/upper bounds for integer variables
% 091019  hkh  New Warm start, input Mx0 - select Mx0 best x_0 values in warm start
% 091019  hkh  Generate all integer combinations. Remove nonfeasible combinations
% 091020  hkh  Handle the case where no x_0 is usable, Inform=9, Exit=999, ExitFlag=4
% 091128  hkh  Handle huge x-range, revise random x_0, limit xInit and use MaxFunc for derivative-free optimization
% 091201  hkh  Set lower tol, use internal derivatives if DFREE=1 for SOL solvers. Given x_0 used 2nd if OK
% 100509  hkh  Safe guard for inf bounds, use x_L,x_U, not x_LL,x_UU if HUGE
% 100511  hkh  Use y, not z, for work array if HUGE=1, avoid length conflicts
% 101126  hkh  Set x_k=min(Prob.x_U,max(Prob.x_L,r.x_k)) x_L,x_U from Prob needed, if |x_k|>1000
% 101126  hkh  Added output from subsolver if PriLevOpt>10. PriLevOpt-10 sent to subsolver
% 101126  hkh  Lower and upper bounds, and constraint values printed if PriLevOpt > 13 (not doc)
% 110127  per  Fix bug that caused one of the integer solutons to be skipped.
% 110411  and  Fix A*empty for linear infeasible problems
% 110723  hkh  Two major loops now run with parfor if Prob.Threads > 1
