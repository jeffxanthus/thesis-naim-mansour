% glcCluster.m
%
% glcCluster solves general constrained mixed integer global optimization
%
% glcCluster is a hybrid algorithm, that is using one of the following DIRECT
% algorithms: glcDirect (default),glcFast or glcSolve, for global search (Step 1). 
% Step 2 is an adaptive clustering algorithm to find a suitable number of clusters, 
% where the best point in each cluster is then used as an initial point for a 
% local search (Step 3).
% The 4th step is to run the DIRECT algoirithm once again, to possible improve.
% If the DIRECT algorithm improves the best point, a local search is 
% finally made as Step 5 with the new best point(s) as starting points.
% For a more detailed algorithm description, see below after USAGE:
%
% glcCluster solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Usage: See below
%
% Calling syntax:
%
% function Result = glcCluster(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)
%
% INPUT PARAMETERS
%
% The following three parameters are either set in Prob.GO, or as extra input
%
% maxFunc1  Number of function evaluations in 1st call to the DIRECT solver.
%           Should be odd number (automatically corrected).
%           Default 200*dim(x) + 1.  May also be set as Prob.GO.maxFunc1
%           Default 500*dim(x) + 1 if any integer variable (black-box MINLP problem)
% maxFunc2  Number of function evaluations in 2nd call to the DIRECT solver.
%           May also be set as Prob.GO.maxFunc2
%           By default maxFunc2 is set to minimum of maxFunc1 and the number
%           of function evaluations left after doing previous phases.
% maxFunc3  If the DIRECT solver is not feasible after maxFunc1 function evaluations,
%           it will be repeatedly called (warm start) doing maxFunc1 function
%           evaluations until maxFunc3 function evaluations reached.
%           May also be set as Prob.GO.maxFunc3
% ProbL     Structure to be used in the local search. By default the same
%           problem structure as in the global search is used, Prob (see below)
%           Using a second structure is important if optimal continuous
%           variables may take values on bounds. the DIRECT solver used for the global
%           search only converges to the simple bounds in the limit, and
%           therefore the simple bounds may be relaxed a bit in the global
%           search. Also, if the global search has difficulty fulfilling
%           equality constraints exactly, the lower and upper bounds may be
%           slightly relaxed. But being exact in the local search.
%           Note that the local search is using derivatives, and can utilize
%           given analytic derivatives. Otherwise the local solver is using
%           numerical derivatives or automatic differentiation.
%           If routines to provide derivatives are given in ProbL,
%           they are used. If only one structure Prob is given, glcCluster
%           uses the derivative routines given in the this structure
%           May also be set as Prob.GO.ProbL
%
% Prob    Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say GLCF
%   FUNCS.c   The routine to compute the nonlinear constraint, say GLCC
%             A call to tomFiles.m or glcAssign.m sets these fields.
%   x_0       Column vector giving an initial point for a local search.
%             If to give a matrix of initial points, use Prob.X0.
%             if [] or just 0s, not used. See Phase 3 in the algorithm
%   X0        Matrix of points (columns), where each is used in a local search
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   PriLevOpt Print Level
%             0 = silent.
%             1 = Some output from each glcCluster phase
%             2 = More output from each phase
%             3 = Adaptive clustering information
%             6 = Use PrintResult( ,1) to print summary from each global and
%                 local run
%             7 = Use PrintResult( ,2) to print summary from each global and
%                 local run
%             8 = Use PrintResult( ,3) to print summary from each global and
%                 local run
%
%   WarmStart If true, >0, glcCluster warm starts the DIRECT solver.
%             The DIRECT solver will utilize all points sampled in last run,
%             from one or two calls, dependent on the success in last run.
%             Note: If WarmStart, the DIRECT solver may not be changed
%             Field Prob.MIP.fIP is set to Result.f_k in WarmDefGLOBAL
%
%   fCut      If initial f(x_0) > fCut, no local optimization is done
%             Default fCut = 1E10
%
%   xEqTol    Tolerance to test if new point x_k already defined as
%             optimum: norm(x_k-xOpt(:,i)) <= xEqTol*max(1,norm(x_k))
%             If test fulfilled x_k is assumed to be too close to xOpt(:,i)
%             Default xEqTol = 1E-5
%
%   fEqTol    Tolerance to test if new f(x_k(1)) is same as (x_k(2))
%             Default fEqTol = 1E-7
%
% ---------------------------------------
% optParam    Structure in Prob, Prob.optParam.
% ---------------------------------------
%             Defines optimization parameters. Fields used:
%  IterPrint  In the global optimization (the DIRECT solver):
%             Print iteration #, # of evaluated points, best f(x) and x for
%             each iteration, where the best point has improved
%             IterPrint = 1 implies PriLevOpt >= 1
%   cTol      Nonlinear constraint tolerance
%   MaxIter   Maximal number of iterations in local search, default 1000
%   MaxFunc   Maximal number of function evaluations, default 10000 (roughly)
%             As 1st step needs MaxFunc = 200*n + 1; then MaxFunc >> 100*n
%   EpsGlob   Global/local weight parameter, default 1E-4
%   fGoal     Goal for function value, if empty not used
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%   eps_x     Convergence tolerance in x. All possible rectangles are
%             less than this tolerance (scaled to (0,1) )
%             See the output field maxTri.
%
% ---------------------------------------
% GO
% ---------------------------------------
%   localSolver Optionally change local solver used ('snopt' or 'npsol' etc.)
%   maxFunc1    See description above
%   maxFunc2    See description above
%   maxFunc3    See description above
%   ProbL       See description above
%
%   DIRECT      DIRECT subsolver, either glcDirect (default),glcFast or glcSolve
%   maxLocalTry Maximal number of local searches from cluster points
%               If <= 0, glcCluster stops after clustering. 
%               Default 30 for continous problems 
%               Default 5 for each unique integer value, if ~isempty(IntVars)
%   minLocalTry Minimal number of local searches from cluster points, Default 1
%   incLocalTry If >0, try to increase number of clusters as much as possible,
%               still being inside [minLocalTry,maxLocalTry], default 0
%               Used for problems which are difficult to find the global minimum for.
%               maxLocalTry might need to be increased as well.
%               Use input Prob.PriLevOpt = 3; to see what clusters are found.
%   maxDistMin  The minimal distance used for clustering, default 0.01
%               See Cluster.maxDist and Cluster.minDist
%   incDist     MaxDist=incDist*MaxDist is used in each adaptive step to decrease 
%               the number of clusters, default 1.1, incDist > 1 
%   decDist     MaxDist=decDist*MaxDist is used in each adaptive step to increase 
%               the number of clusters, default 0.9, decDist < 1
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP.
% ---------------------------------------
%             Defines mixed-integer optimization parameters. Fields used:
%   IntVars:  
%             If empty, all variables are assumed non-integer 
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
%             It is advised to number the integer values as the first
%             variables, before the continuous. The tree search will then
%             be done more efficiently.
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  c_k      Nonlinear constraints values at x_k
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ConstrEv Number of constraint function evaluations
%  maxTri   Maximum size of any triangle
%  ExitText Text string giving ExitFlag and Inform information
%
%  -------------------------------------------------------
%  Result.Cluster, subfield with clustering information
%  -------------------------------------------------------
%  Cluster.x_k      Matrix with best cluster points
%  Cluster.f_k      Row vector with f(x) values for each column in Cluster.x_k
%  Cluster.maxDist  maxDist used for clustering. Points with the distance to
%                   nearest neighbor (minDist) closer than maxDist are
%                   clustered
%  Cluster.minDist  Vector with the minimal distance from each point to any
%                   neighbor
%  Cluster.minIdx   Vector with indices for the points used in clustering
%
% USAGE:
% The algorithm has the following five phases:
% Phase 0: Do local search starting with the column vector given in Prob.x_0 
%          (if not [] or all 0s) and each of the points given as columns in 
%          Prob.X0 (if defined). Also start from adjusted center point.
%          Best point (xMin,fMin) is input to DIRECT routine in Phase 1.
% Phase 1: Run DIRECT solver maxFunc1 function value trials
%          If DIRECT never finds a feasible point, it will run warm starts
%          with maxFunc1 function evaluations until totally maxFunc3 evaluations
%          If DIRECT never finds a feasible point, and maxFunc is reached
%          then glcCLuster will stop
% Phase 2: Apply a clustering algorithm on all sampled points by DIRECT
%          The algorithm finds a set of nPnt clusters. The point with the
%          lowest function value in each cluster is selected.
%          If the number of clusters are below minLocalTry or above maxLocalTry,
%          an iteration to find a suitable number of clusters are done by 
%          changing the distance parameter maxDist, see incDist and decDist
%          If incLocalTry > 0, an iteration decreasing maxDist (and increasing
%          the number of clusters) is tried and the maximal number of 
%          clusters <= maxLocalTry is used 
%          When doing the clustering, all point with distance >= mean(all distances)
%          are considered far away and not used. The lowest f(x) of these, Fout, 
%          is found. If Fout < any of the best cluster points, this (these) points
%          are added to the set of initial points.
% Phase 3: Do local search with of the nPnt best cluster points as
%          initial starting value. Most likely the local search will find
%          the global optimum, if there are not too many local minima
% Phase 4: If the best point in the local searches in Phase 3 is better 
%          than the best point found in the global search in Phase 1, this 
%          new point (xMin,fMin) is added as input to the DIRECT solver and a warm
%          start run with DIRECT doing maxFunc2 function trials is done.
% Phase 5: Apply clustering algorithm as in Phase 2, generating nPnt2 clusters, nPnt2 >= nPnt.
%          Select points not too close to previous initial points.
% Phase 6: Local search from each of the nPnt2 points with the best function value
%          If the local search improves the best point found, then this point (xMin,fMin)
%          could be used as input in Prob.xIP = xMin, Prob.fIP = fMin if the user does a further
%          warm start of the DIRECT solver. 
%          A warm start of glcCluster is also possible.
%
% Let the name of the problem be "GLCF Test"
% The function GLCF is best written as
%     function f = GLCF(x, Prob)
% Then any information, say u and W is easily sent to GLCF (and GLCC) using
% the Prob structure. See the example below.
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%    Prob   = glcAssign('GLCF',x_L,x_U,'GLCF Test',A,b_L,b_U,'GLCF',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown below
%
%    Prob.GO.maxFunc1 = 1000; % Number of initial function evals to DIRECT solver
%    Prob.GO.maxFunc2 = 200;  % Number of final   function evals to DIRECT solver
%    Prob.GO.maxFunc3 = 2000; % Total number of calls in Phase 1, if infeasible
%    Prob.GO.ProbL = ProbL;   % Set the structure for the local search
%
%    Prob.optParam.MaxFunc = 10000; % Must be sufficiently large for all local
%                                   % searches + maxFunc1 + maxFunc2  and
%                                   % if nonfeasible also: + (maxFunc3-maxFunc1)
%    Prob.optParam.MaxIter = 1000; % Max number of iterations in each local
%                             % search, normally much less are used
%    Prob.PriLevOpt = 1;      % Some printing from each algorithm phase
%
%    If there are integer variables, they may be set as additional input
%    to glcAssign, or directly as the next line shows:
%    Prob.MIP.IntVars = [1 3];  % 1st and third variables are integers
%
% A second structure ProbL for the local search is best setup using
% conAssign.
%
% Driver call, including printing with level 2:
%    Result = tomRun('glcCluster',Prob,2);
%
% Direct solver call:
%    Result = glcCluster(Prob);
% or with all input parameters as
%    Result = glcCluster(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)
% Then print the results with call
%    PrintResult(Result);
%
% The user function GLCF is written as
%
%    function f = GLCF(x, Prob)
%    u = Prob.user.u; W = Prob.user.W;
%    f = "some function of x, u and W"
%
% It is also possible to use the function format
%    function f = GLCF(x)
% but then any additional parameters must be sent as global variables.
%
% The user function GLCC, computing the nonlinear constraints, is written as
%
%    function c = GLCC(x, Prob)
%    u = Prob.user.u; W = Prob.user.W;
%    c = "some vector function of x, V and W"
%
% Note! If GLCF has the 2nd input argument Prob, also GLCC must have that.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 3, 2001.    Last modified July 28, 2011.

function Result = glcCluster(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)

if nargin < 5
   ProbL = [];
   if nargin < 4
      maxFunc3 = [];
      if nargin < 3
         maxFunc2 = [];
         if nargin < 2
            maxFunc1 = [];
            if nargin < 1
               error('glcCluster needs input structure Prob');
            end
         end
      end
   end
end

global NARG NLP_f NLP_c NLP_g NLP_dc NLP_r NLP_J NLP_x

solvType=checkType('glc');

if isfield(Prob.GO,'DIRECT')
   DIRECT = Prob.GO.DIRECT;
else
   DIRECT = [];
end
if strcmpi(DIRECT,'glcSolve')
   DIRECT     = 'glcSolve';
   DIRECTFILE = 'glcSave';
   OLD        = 1;
elseif strcmpi(DIRECT,'glcFast')
   DIRECT     = 'glcFast';
   DIRECTFILE = 'glcFastSave';
   OLD        = 1;
else
   DIRECT     = 'glcDirect';
   DIRECTFILE = 'glcDirectSave';
   OLD        = 0;
end

Prob=ProbCheck(Prob,DIRECT,solvType);

if isempty(Prob.x_L) | isempty(Prob.x_U)
   disp('glcCluster requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'glcCluster requires both lower and upper variable bounds';
   Result=endSolve(Prob,Result);
   return;
end

Prob = iniSolve(Prob,solvType,0,0);

mFunc = Prob.optParam.MaxFunc;       % Maximum number of function evaluations
mIter = Prob.optParam.MaxIter;       % Number of iterations
% Safeguard
if mIter < 0
   mIter = 10000;
end
if mFunc < 0
   mFunc = 10000;
end

PriLev    = Prob.PriLevOpt;          % Print level
if PriLev > 4
   PriLevCl = 2;
elseif PriLev <= 0
   PriLevCl = 0;
else
   PriLevCl = 1;
end
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
xEqTol    = DefPar(Prob,'xEqTol',1E-5);
fEqTol    = DefPar(Prob,'fEqTol',1E-7);
fCut      = DefPar(Prob,'fCut',1E10);
cTol      = Prob.optParam.cTol;      % Tolerance in constraints
if IterPrint
   PriLev = max(1,PriLev);
end

if isinf(fGoal), fGoal = -1E300; end


if isfield(Prob.GO,'localSolver')
   localSolver   = deblank(Prob.GO.localSolver);
else
   localSolver = GetSolver('con',0,0);
end

if isfield(Prob.GO,'minLocalTry')
   minLocalTry   = Prob.GO.minLocalTry;
else
   minLocalTry = [];
end
if isempty(minLocalTry), minLocalTry = 1; end

if isfield(Prob.GO,'incLocalTry')
   incLocalTry   = Prob.GO.incLocalTry;
else
   incLocalTry = [];
end
if isempty(incLocalTry), incLocalTry = 1; end

if isfield(Prob.GO,'incDist')
   incDist   = max(0,Prob.GO.incDist);
else
   incDist = [];
end
if isempty(incDist), incDist = 1.1; end
if incDist <= 1, incDist = 1.1; end

if isfield(Prob.GO,'decDist')
   decDist   = max(0,Prob.GO.decDist);
else
   decDist = [];
end
if isempty(decDist), decDist = 1.1; end
if decDist >= 1, decDist = 0.9; end


Result                 = ResultDef(Prob);
Result.Solver          = 'glcCluster';
Result.SolverAlgorithm = ...
 ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];

% Do a preSolve, hopefully shrinking the box
Prob.PriLevOpt = Prob.PriLevOpt-1;
% Decide if to do a presolve or not, now set to do a presolve
if 1
   Prob2 = preSolve(Prob); 
   Prob.PriLevOpt = Prob.PriLevOpt+1;
   % Pick up input parameters from the Prob structure:
   x_L = Prob2.x_L(:);   % Lower bounds
   x_U = Prob2.x_U(:);   % Upper bounds
else
   x_L = Prob.x_L(:);   % Lower bounds
   x_U = Prob.x_U(:);   % Upper bounds
end

% Check for Inf and set to lower values.
x_L(isinf(x_L)) = -10000;
x_U(isinf(x_U)) =  10000;

x_0      = Prob.x_0;
X0       = DefPar(Prob,'X0',[]);
FX0      = [];
XX       = [];
Prob.x_L = x_L;
Prob.x_U = x_U;
x_D      = x_U - x_L;
ixD      = find(x_D > 0);
n        = Prob.N;

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('glcCluster: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars      = find(IV);
nI           = length(IntVars);
Reals        = find(~IV);
nR           = length(Reals);
% WarmStart  = DefPar(Prob,'WarmStart',0);


if length(IntVars) == n
   fprintf('glcCluster: Do not call glcCluster for PURE integer programs\n');
   fprintf('            ');
   fprintf('No point in clustering and local search for pure IP\n');
   fprintf('            ');
   fprintf('Call glcDirect, glcSolve, or some other IP solver instead\n');
   error('glcCluster: Wrong type of optimization problem');
end
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

fCut           = DefPar(Prob,'fCut',1E10);

if isempty(Prob.FUNCS.g)
   DFREE       = 1;
else
   DFREE       = 0;
end

if DFREE
   if SOLSNOPT | SOLNPSOL
      Prob.SOL.optPar(9)  = 1E-4;
      Prob.SOL.optPar(10) = 1E-4;
      Prob.SOL.optPar(40) = 0;
      Prob.NumDiff        = 6;
      if Prob.mNonLin > 0
         Prob.ConsDiff    = 6;
      end
      %%%Prob.SOL.optPar(22) = 0.1;
      % Comment following three lines for now, should check more best choices
      %Prob.SOL.PrintFile  = 'SOL.txt';
      %Prob.SOL.optPar(1)  = 111111;
      %Prob.SOL.optPar(5)  = 1;
      %%%Prob.SOL.optPar(39) = 0;
   end
end

if isfield(Prob.GO,'maxLocalTry')
   maxLocalTry   = Prob.GO.maxLocalTry;
else
   maxLocalTry = [];
end
if isempty(maxLocalTry) 
   if isempty(IntVars)

      if DFREE
         if mFunc < 5000
            maxLocalTry = 5; 
         elseif mFunc < 10000
            maxLocalTry = 10; 
         else
            maxLocalTry = 20; 
         end
      else
         if mFunc < 5000
            maxLocalTry = 10; 
         elseif mFunc < 10000
            maxLocalTry = 20; 
         else
            maxLocalTry = 30; 
         end
      end
   else
      maxLocalTry = 5; 
   end
end
if isfield(Prob.GO,'maxDistMin')
   maxDistMin   = Prob.GO.maxDistMin;
else
   maxDistMin = [];
end
if isempty(maxDistMin) 
   if isempty(IntVars)
      maxDistMin = 0.01; 
   else
      maxDistMin = 0.001; 
   end
end

if isempty(maxFunc1)
   if isfield(Prob.GO,'maxFunc1')
      maxFunc1 = Prob.GO.maxFunc1;
   end
end
if isempty(maxFunc1) 
   if isempty(IntVars)
      if DFREE
         maxFunc1 = min(0.2*mFunc,max(200*n,0.1*mFunc));
      else
         maxFunc1 = min(0.3*mFunc,max(200*n,0.2*mFunc));
         % HKH Should test for best maxFunc1 choice
         %maxFunc1 = 200*n + 1;    % Odd number wanted :)
      end
      % Odd number wanted :)
      if mod(maxFunc1,2) == 0, maxFunc1 = maxFunc1+1; end
   else
      maxFunc1 = 200*n + 1;    % Odd number wanted :)
   end
else
   % Make sure it is odd
   if mod(maxFunc1,2) == 0, maxFunc1 = maxFunc1+1; end
end

if isempty(maxFunc2)
   if isfield(Prob.GO,'maxFunc2')
      maxFunc2 = Prob.GO.maxFunc2;
   end
end
if isfield(Prob.GO,'maxFunc3')
   maxFunc3 = max(maxFunc1,Prob.GO.maxFunc3);
elseif ~isempty(maxFunc3)
   maxFunc3 = max(maxFunc1,maxFunc3);
else
   maxFunc3 = min(mFunc,3*maxFunc1);
end

maxFunc = maxFunc1;
maxIter = maxFunc;

if (maxFunc >= mFunc)
   fprintf('maxFunc1 == MaxFunc 1st step = %d\n',maxFunc);
   fprintf('That is not possible, because total max is %d!\n',mFunc);
   fprintf(' ');
   error('glcCluster - Illegal input!!!')
end
totfEval = 0;
totcEval = 0;
gradEv   = 0;
HessEv   = 0;
useProb  = 0;
if isempty(ProbL)
   if isfield(Prob.GO,'ProbL')
      ProbL = Prob.GO.ProbL;
      %ProbL = ProbCheck(ProbL,localSolver);
   end
   if isempty(ProbL)
      useProb = 1;
      ProbL = Prob;
   end
else
   %ProbL = ProbCheck(ProbL,localSolver);
end
Prob.solvType = solvType;
% Set print levels
tomPrint        = PriLev - 5;
Prob.PriLevOpt  = 0;
ProbL.PriLevOpt = 0;
ProbL.WarmStart = 0;
if ~useProb
   N =ProbL.N;
   M = Prob.mLin + Prob.mNonLin;
   if isfield(ProbL,'optParam')
      if isempty(ProbL.optParam) | ProbL.probType ~= 3
         ProbL.optParam=optParamDef(localSolver,3,N,N,M);
      else
         ProbL.optParam=optParamSet(ProbL.optParam,localSolver,3,N,N,M);
      end
   else
      ProbL.optParam=optParamDef(localSolver,3,N,N,M);
   end
else
   probType       = 3;
   N              = ProbL.N;
   M              = Prob.mLin + Prob.mNonLin;
   ProbL.optParam = optParamDef(localSolver,probType,N,N,M);
end
% optP=ProbL.optParam
% Check if bounds huge, or squeezed
ix = x_D(Reals) >= 1000; 
if any(ix) 
   HUGE = 1;
else
   HUGE = 0;
end

cDev  = [];
if Prob.WarmStart == 1
   load(DIRECTFILE,'F');
   fPrev = length(F);
   if isfield(Prob.MIP,'fIP')
      fMin  = Prob.MIP.fIP;
      xMin  = Prob.MIP.xIP;
   else
      fMin  = inf;
   end
   % fMin  = min(F);
   % nX0              = 0;
   % FX0              = [];
   % XX               = ones(n,0);
   % Result.Cluster.maxDist = maxDistV;
   % Result.Cluster.minDist = minlV;
   % Result.Cluster.minIdx  = v2V;
   % Result.Cluster.x_k     = x_k;
   % Result.Cluster.f_k     = Fv;
   % Result.Cluster.c_k     = c_k;
   % Saved local tries
   % Result.Cluster.FX0     = FX0;
   % Result.Cluster.X0      = X0;
   % Result.Cluster.XX      = XX;
   % Pick up the points evaluated previous ru
   FX0    =  Prob.Cluster.FX0;
   X0     =  Prob.Cluster.X0;
   XX     =  Prob.Cluster.XX;
   nX0    = size(XX,2);
elseif maxLocalTry > 0
   % ==============================================
   % Phase 0 Local search from x_0, X0 and xS
   % ==============================================
   fMin  = Inf;
   xMin  = [];
   fPrev = 0;
   NARG  = []; NLP_x = []; NLP_f = []; NLP_c = []; NLP_g = []; NLP_dc = []; NLP_r = []; NLP_J = []; 
   if ~isempty(X0)
      if PriLev > 0
         fprintf('   %d points in Prob.X0.',size(X0,2))
      end
   end
   if all(x_0 == 0), x_0 = []; end
   if ~isempty(x_0)
      X0 = [x_0,X0];
      if PriLev > 0
         fprintf('   1 point in Prob.x_0.')
      end
   end
   % Find base point
   if HUGE
      % Find a base point close to 0
      ixH           = find(ix);
      ixI           = find(~ix);
      xB            = zeros(nR,1);
      if nI > 0
         xLL        = x_L(Reals);
         xUU        = x_U(Reals);
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
      z(ixI)        = (rand(length(ixI),1)-0.5).*xDD(ixI);
      Fac           = 10*Fac;
      %if Fac > 1, Fac = 10^(-Pow); end
      z(ixH)        = Fac*randn(length(ixH),1).*xDD(ixH);
      z(ixHp)       = abs(z(ixHp));
      z(ixHm)       = -abs(z(ixHm));
      %X0(RealVars,P)  = xB + z;
      xB0           = x_L;
      xB0(Reals)    = xB + z;
   else
      xB0           = x_L;
      xB0(Reals)    = (x_L(Reals)+x_U(Reals))/2;
   end
   X0               = [X0,xB0];
   if PriLev > 0
      fprintf('\n')
   end
   nX0              = size(X0,2);
   FX0              = zeros(2,nX0);
   XX               = Inf*ones(n,nX0);

   for i = 1:nX0
       FX0(1,i)  = nlp_f(X0(:,i), Prob);         % Compute f(x_0) for each x_0
       totfEval  = totfEval + 1;
   end 
   [F00,ix]      = sort(FX0(1,:));
   ixBest        = 0;
   for j=1:nX0
     i           = ix(j);
     NLP_x       = X0(:,i);
     NLP_f       = FX0(1,i);
     ProbL.x_0   = NLP_x;
     ProbL.f_0   = NLP_f;
     if NLP_f >= fCut
        FX0(2,i) = inf;
     else
       if ~isempty(IntVars)
          x00    = max(x_L(IntVars),min(x_U(IntVars),ProbL.x_0(IntVars)));
          ProbL.x_L(IntVars) = x00;
          ProbL.x_U(IntVars) = x00;
       end
       rl          = tomRunFast(localSolver,ProbL);
       PrintResult(rl,tomPrint);
       totfEval    = totfEval + rl.FuncEv;
       totcEval    = totcEval + rl.ConstrEv;
       gradEv      = gradEv + rl.GradEv;
       HessEv      = HessEv + rl.HessEv;
       if M > 0
          % Check if point feasible enough to be accepted
          h_L1     = consviolation(rl, 0);
          if h_L1 <= M*ProbL.optParam.cTol
             f_k      = rl.f_k;
             if PriLev > 2
                fprintf('   ')
                fprintf('Local point %3d accepted, f(x) = %20.10f.',i,f_k);
                fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
                fprintf(' sum(|constr|)%18.10f',h_L1);
                if isempty(IntVars)
                   fprintf('\n');
                else
                   xprinti(rl.x_k(IntVars),' IV:');
                end
             end
          else
             f_k = Inf;
             if PriLev > 3
                fprintf('   ')
                fprintf('Local point %3d NOT OK,   f(x) = %20.10f.',i,rl.f_k);
                fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
                fprintf(' sum(|constr|)%18.10f',h_L1);
                if isempty(IntVars)
                   fprintf('\n');
                else
                   xprinti(rl.x_k(IntVars),' IV:');
                end
             end
          end
       else
          f_k = rl.f_k;
          if PriLev > 2
             fprintf('   ')
             fprintf('Local point %3d, f(x) = %20.10f.',i,f_k);
             fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
             if isempty(IntVars)
                fprintf('\n');
             else
                xprinti(rl.x_k(IntVars),' IV:');
             end
          end
       end
       FX0(2,i) = f_k;
       XX(:,i)  = rl.x_k;
       if f_k <= fMin
          if f_k == fMin
             xMin   = [xMin,rl.x_k];
             cDev   = [cDev,rl.c_k];
             ixBest = [ixBest,i];
          else
             xMin   = rl.x_k;
             cDev   = rl.c_k;
             ixBest = i;
          end
          fMin      = f_k;
          CMin      = ProbL.x_L;
          CMin(ixD) = (rl.x_k(ixD) - x_L(ixD))./x_D(ixD);
          if ~isempty(IntVars)
             IMin   = ProbL.x_L(IntVars);
          end
          if PriLev > 1 
             fprintf('   ')
             fprintf('Step %3d. Pnt %3d. ',j,i)
             fprintf(' New Best %18.16f',f_k)
             if isempty(IntVars)
                fprintf('\n');
             else
                xprinti(rl.x_k(IntVars),'. IV:');
             end
          end
       end
       if ~isempty(IntVars)
          ProbL.x_L(IntVars) = x_L(IntVars);
          ProbL.x_U(IntVars) = x_U(IntVars);
       end
       if PriLev > 5 
          fprintf('   ')
          fprintf('End of local search # %d out of %d\n',i,size(X0,2)); 
          fprintf('                      ')
          fprintf('----------------------------\n\n'); 
       end
       if totfEval > mFunc, break; end
     end
   end

   if ixBest(1) > 0
      if PriLev > 0
         fprintf('glcCluster Ph-0: ')
         fprintf('Local solver %s. ',localSolver)
         fprintf('Local try')
         fprintf(' %d',ixBest)
         fprintf('/%d: ',nX0)
         fprintf('fMin %18.16f.',fMin)
         if isempty(IntVars)
            fprintf('\n');
         else
            xprinti(xMin(IntVars,1),' IV:');
            for j=2:size(xMin,2)
                xprinti(xMin(IntVars,j),' IV:');
            end
         end
      end
   else
      fMin         = Inf;
      if PriLev > 0
         fprintf('glcCluster Ph-0: ')
         fprintf('No feasible point found in %d local tries.\n',nX0)
      end
   end
end

% ==========================================
% Phase 1 First DIRECT run
% ==========================================
Prob.optParam.MaxFunc   = maxFunc;
Prob.optParam.MaxIter   = maxFunc;   % MaxIter is of less interest

% Set best point in Prob for input to DIRECT solver
if ~isinf(fMin)
   Prob.MIP.fIP         = fMin;
   Prob.MIP.xIP         = xMin(:,end);
end

% 1st step: Run glc with maxFunc evaluations
Prob.optParam.IterPrint = PriLev > 5;
r                       = tomRunFast(DIRECT,Prob);
PrintResult(r,tomPrint);
Prob.optParam.IterPrint = 0;

load(DIRECTFILE,'feasible'); 

totfEval                = totfEval + r.FuncEv;
totcEval                = totcEval + r.ConstrEv;

% Avoid more testing of the structure
Prob.MENU = 1;

if ~feasible
   Prob.optParam.MaxFunc = maxFunc-1;
   ExitFlag = r.ExitFlag;
   while ~feasible & (totfEval < maxFunc3) & (ExitFlag == 0 | ExitFlag == 7)
      % Must try to get feasible
      Prob.WarmStart = 1;
      r = tomRunFast(DIRECT,Prob);
      PrintResult(r,tomPrint);
      totfEval = totfEval + r.FuncEv;
      totcEval = totcEval + r.ConstrEv;
      ExitFlag = r.ExitFlag;
      load(DIRECTFILE,'feasible'); 
   end
   Prob.WarmStart = 0;  % Reset WarmStart
end
if ~feasible & totfEval >= mFunc
   Result = r;
   Result.Solver          = 'glcCluster';
   Result.SolverAlgorithm = ...
    ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of major iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 7;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev) ...
                      ' f(x) > max f(x) evaluations. No feasible point found'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('Max f(x) evaluations reached. No feasible point found\n'); end
   return
end

Prob.WarmStart = 0;  % No WarmStart in the local solvers

% C and F could have been returned directly from glcDirect/glcSolve

if OLD
   load(DIRECTFILE,'C','glcfMin','F','ignoreIdx','fMinIdx','fMinEQ','G'); 
else
   load(DIRECTFILE,'C','glcfMin','F','ign','fMinIdx','fMinEQ','G'); 
   ignoreIdx = ign;
end
if PriLev > 0
   fprintf('glcCluster Ph-1: ')
   fprintf('Tried %4d. ',totfEval)
   if feasible
      fprintf('Best feasible global point glcfMin:')
      fprintf('%20.10f.',glcfMin)
      if fMinEQ > 0
         fprintf('                      ')
         fprintf('f(x)-fMinEQ = %14.10f ',glcfMin-fMinEQ);
         fprintf('Eq.Cons Dev fMinEQ %12.10f\n',fMinEQ);
      else
         fprintf('\n');
      end
   else
      fprintf('No feasible point found! ')
      if fMinEQ > 0
         fprintf('Eq.Cons Dev fMinEQ %12.10f\n',fMinEQ);
      else
         fprintf('\n');
      end
   end
end
% Update fMin if glcfMin from DIRECT better than Phase 0 solution
if glcfMin < fMin & feasible
   fMin = glcfMin;
   xMin = r.x_k;
   cDev = r.c_k;
end
% Use fMinEQ nonzero to flag for local search to accept any feasible point
if ~feasible & fMinEQ == 0 & isinf(fMin)
   fMinEQ = 1000000;
end
% ==========================================
% Phase 2 First clustering
% ==========================================

% HKH New code to handle integer variables and feasibility

if M == 0
   iPnt  = ones(1,length(F));
   ix    = [1:length(F)];
else
   iPnt  = double(all(G <= cTol));
   ix    = find(iPnt > 0);
   if isempty(ix)
      cTol = 100*cTol;
      if PriLev > 0
         fprintf('   Increasing cTol %f\n',cTol);
      end
      iPnt  = double(all(G <= cTol));
      ix    = find(iPnt > 0);
   end
   if isempty(ix)
      % Desperate - avoid crash
      iPnt  = ones(1,length(F));
      ix    = [1:length(F)];
   end
   % Must set nonused elements to <0, because integer keys from 0,1,...
   iPnt(iPnt==0) = -inf;
end

F0        = [];
n0        = 0;
c_0       = [];
nOK       = 0;
minlV0    = [];
v2V0      = [];
maxDistV0 = [];

if isempty(IntVars)
   iVal     = 1;
   nVal     = 1;
   Fout     = Inf*ones(nVal,1);
   Cout     = zeros(n,nVal);
else
   xD       = 1 + x_D(IntVars);
   iPnt(ix) = C(IntVars(1),ix);
   nn       = 1;
   for i = 2:length(IntVars)
       nn       = nn*xD(i-1);
       iPnt(ix) = iPnt(ix) + nn*C(IntVars(i),ix);
   end
   iVal     = unique(iPnt(ix));
   nVal     = length(iVal);

   Fout     = Inf*ones(nVal,1);
   Cout     = zeros(n,nVal);
   Elements = prod(xD);
   if Elements > nVal
      % Check the integer combinations not in iVal
      %iNot  = double(~all(G <= cTol));
      if ~isempty(Prob.A)
         AI = Prob.A(:,IntVars);
         if all(all(AI==0))
            AI = [];
         else
            Ap = max(0,Prob.A(:,Reals));
            Am = min(0,Prob.A(:,Reals));
         end
      else
         AI    = [];
      end
      iNot     = double(isinf(iPnt));
      ix       = find(iNot > 0);
      iNot(ix) = C(IntVars(1),ix);
      nn       = 1;
      for i = 2:length(IntVars)
          nn       = nn*xD(i-1);
          iNot(ix) = iNot(ix) + nn*C(IntVars(i),ix);
      end
      iVal2     = unique(iNot(ix));
      nVal2     = length(iVal2);
      maxDistV0 = NaN*ones(nVal2,1);
      %nPnts    = 2*(Elements-nVal); % Too big number
      nPnts     = max(2,min(maxLocalTry*nVal2,mFunc));
      F0        = zeros(nPnts,1);
      c_0       = zeros(n,nPnts);
      for i = 1:nVal2
          j = iVal2(i);
          if ~(any(j==iVal))
             ix    = find(iNot == j);
             ixLen = length(ix);
             if isempty(AI)
                isOK = 1;
             else
                % Check if this integer values is possibly feasible
                x_I = round(tomsol(9,x_L(IntVars),C(IntVars,ix(1)),x_D(IntVars)));   
                AIxI = AI*x_I;
                rhs = Prob.b_U-Ap*x_L(Reals)-Am*x_U(Reals);
                lhs = Prob.b_L-Ap*x_U(Reals)-Am*x_L(Reals);
                isOK = all(AIxI <= rhs+1E-6) & all(AIxI >= lhs-1E-6);
                %if ~isOK, keyboard, end
                if PriLev > 0
                   fprintf('   isOK %d: ',isOK);
                   xprinti(x_I,'IV: ');
                end
             end
             if isOK
                nOK                  = nOK + 1;
                % Call cluster algorithm
                [FvNew,cNew,Fout(nOK),Cout(Reals,nOK),maxDist,minl,v2] = clustAlg(...
                   C(Reals,ix),F(ix),PriLevCl,minLocalTry,maxLocalTry,incLocalTry,...
                   incDist,decDist,maxDistMin);

                nPnt                    = length(FvNew);
                F0(n0+1:n0+nPnt)        = FvNew;
                if ~isempty(IntVars)
                   c_0(IntVars,n0+1:n0+nPnt) = C(IntVars,ix(1:nPnt));
                   Cout(IntVars,nOK)         = C(IntVars,ix(1));
                end
                c_0(Reals,n0+1:n0+nPnt) = cNew;
                n0                      = n0 + nPnt;
                maxDistV0(i)            = maxDist;
                if isempty(minlV0)
                   minlV0       = minl;
                else
                   minlV0       = [minlV0,0,minl];
                end
                if isempty(v2V0)
                   v2V0         = v2;
                else
                   v2V0         = [v2V0,0,v2];
                end
                %if size(G,1) == 1
                %   n0  = n0+1;
                %   z = F(ix);
                %   [F0(n0) Iout] = min(z(:)'+max(0,G(:,ix)));
                %   F0(n0)        = F(ix(Iout));
                %   c_0(:,n0)     = C(:,ix(Iout));
                %elseif size(G,1) > 1
                %   n0  = n0+1;
                %   z = F(ix);
                %   [F0(n0) Iout] = min(z(:)'+max(0,sum(G(:,ix))));
                %   F0(n0)        = F(ix(Iout));
                %   c_0(:,n0)     = C(:,ix(Iout));
                %end
                %n0  = n0+1;
                %[F0(n0) Iout] = min(F(ix));
                %c_0(:,n0)     = C(:,ix(Iout));
             end
             %end
          end
      end
      if n0 < nPnts
         F0   = F0(1:n0);
         c_0  = c_0(:,1:n0);
      end
      if PriLev > 5
         xprint(F0,'F0: ')
         xprint(c_0,'c_0:',[],n)
      end
   end

   if PriLev > 0
      fprintf('   ')
      fprintf('Integer combinations %d.',Elements);
      fprintf(' Feasible: %d.',nVal);
      fprintf(' InFeasible:  %d.',nOK);
      fprintf(' Integer inFeasible: %d.',Elements-nVal-nOK);
      fprintf('\n');
   end
end

% Update points to ignore, not considered for trisection 
% Later all points that are clustered are set to be ignored
ignore            = zeros(length(F),1);
ignore(ignoreIdx) = 1;
% Best single point sampled, not included in any cluster
nPnts             = min(1000,maxLocalTry*nVal);
Fv                = zeros(nPnts,1);
c_k               = zeros(n,nPnts);
nv                = 0;
maxDistV          = NaN*ones(nVal,1);
maxDist           = NaN;
minlV             = [];
v2V               = [];
i0                = nOK;

for i = 1:nVal
   % Indices with fixed integer combination corresponding to key iVal(i)
   iKey = iVal(i);
   ix   = find(iPnt == iKey);
   % Call cluster algorithm
   [FvNew,cNew,Fout(i0+i),Cout(Reals,i0+i),maxDist,minl,v2] = clustAlg(...
      C(Reals,ix),F(ix),PriLevCl,minLocalTry,maxLocalTry,incLocalTry,...
      incDist,decDist,maxDistMin);

   nPnt                    = length(FvNew);
   Fv(nv+1:nv+nPnt)        = FvNew;
   if ~isempty(IntVars)
      c_k(IntVars,nv+1:nv+nPnt) = C(IntVars,ix(1:nPnt));
   end
   c_k(Reals,nv+1:nv+nPnt) = cNew;
   nv                      = nv + nPnt;
   maxDistV(i)             = maxDist;
   if isempty(minlV)
      minlV         = minl;
   else
      minlV         = [minlV,0,minl];
   end
   if isempty(v2V)
      v2V         = v2;
   else
      v2V         = [v2V,0,v2];
   end
end

%if nVal > 0
%   if isempty(IntVars)
%      % ix = find(Fout <= max(Fv(1:nv)));
%      ix = [1:nOK+nVal];
%   else
%      ix = [1:nOK+nVal];
%   end
%else
%   ix = [1:nOK];
%end

% Safer to always use the best far away point
ix = [1:nOK+nVal];
nx = length(ix);
if nx > 0
   % Add start from best non-clustered points
   Fv(nv+1:nv+nx)    = Fout(ix);
   c_k(:,nv+1:nv+nx) = Cout(:,ix);
end
nPnt = nv + nx;
if n0 > 0
   % Add start from other infeasible integer combinations
   Fv(nPnt+1:nPnt+n0)    = F0;
   c_k(:,nPnt+1:nPnt+n0) = c_0;
   nPnt                  = nPnt + n0;
end
if nPnt < nPnts
   % Shrink to used size
   Fv  = Fv(1:nPnt);
   c_k = c_k(:,1:nPnt);
end
if PriLev > 4
   xprint(Fv,'Fv: ')
   xprint(c_k,'c_k:',[],n)
end
nTot         = nPnt;
[c_k I J]    = unique(c_k','rows');
c_k          = c_k';
Fv           = Fv(I); 
nPnt         = size(c_k,2);
x_k          = tomsol(9,x_L,c_k,x_D);
x_k(IntVars) = round(x_k(IntVars));

if PriLev > 0
   fprintf('   ')
   fprintf('Total pnts  %d. ',nTot)
   fprintf('Cluster pnts %d. ',nv)
   fprintf('Single  pnts %d. ',nVal)
   if ~isempty(IntVars)
      fprintf('Infeasible cluster pnts %d.',n0)
      fprintf('Single nonfeasible pnts %d. ',nOK)
   end
   if nPnt ~= nTot
      fprintf(' Unique pnts %d!',nPnt)
   end
   fprintf('\n')
end

if isempty(maxFunc2)
   left = mFunc - totfEval;
else
   left = mFunc - totfEval - maxFunc2;
end

if left < nPnt
   % Too few function values left, user must increase MaxFunc
   Result = r;
   Result.Solver          = 'glcCluster';
   Result.SolverAlgorithm = ...
    ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
   Result.Cluster.maxDist = maxDistV;
   Result.Cluster.minDist = minlV;
   Result.Cluster.minIdx  = v2V;
   Result.Cluster.x_k     = x_k;
   Result.Cluster.f_k     = Fv;
   Result.Cluster.c_k     = c_k;
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 1;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev)  ...
                      ' f(x). Max f(x) evaluations reached.'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('   Max f(x) evaluations reached\n'); end
   return
end
if maxLocalTry <= 0
   % Only clustering should be made
   Result = r;
   Result.Solver          = 'glcCluster';
   Result.SolverAlgorithm = ['Constrained DIRECT with ' DIRECT '/ Cluster / NO Local Search'];
   Result.Cluster.maxDist = maxDistV;
   Result.Cluster.minDist = minlV;
   Result.Cluster.minIdx  = v2V;
   Result.Cluster.x_k     = x_k;
   Result.Cluster.f_k     = Fv;
   Result.Cluster.c_k     = c_k;
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 0;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev)  ...
                      ' f(x). Found ' num2str(nPnt) ' cluster points'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('   Clustering done, %d points found\n',nPnt); end
   return
end

% ==============================================
% Phase 3 Local search from best cluster points
% ==============================================
%Prob.optParam.MaxIter=left/nPnt;
Prob.optParam.MaxIter=mIter;
Prob.optParam.MaxFunc=left/nPnt;

X    = [];
Floc = [];

nxk = size(x_k,2);
if PriLev > 0
   fprintf('glcCluster Ph-2: ')
   fprintf('%d initial points from clustering.',nxk)
   fprintf('\n')
end
Fxk = zeros(2,nxk);
XK  = zeros(n,nxk);
if 0 & fMinEQ > 0
   for i = 1:nxk
       Fxk(1,i) = nlp_f(x_k(:,i), ProbL);         % Compute f(x_0) for each x_0
       totfEval = totfEval + 1;
   end 
else
   Fxk(1,:) = Fv';
end


% Remove same/similar initial values
remove = 0;
for P = 1:nX0
    f0   = FX0(1,P);
    if ~isnan(f0)
       %ix   = find(abs(f0-Fxk(1,ix)) < fEqTol*max(1,abs(f0)));
       ix   = find(abs(f0-Fxk(1,:)) < fEqTol*max(1,abs(f0)));
       if ~isempty(ix)
          iE              = checkEq(X0(:,P),x_k(:,ix),xEqTol);
          if ~isempty(iE)
             remove       = remove + length(iE);
             % HKH New
             %F0(1,ix(iE)) = NaN;
             Fxk(1,ix(iE)) = NaN;
          end
       end
    end
end
Fxk(1,Fxk(1,:) >= fCut) = NaN;                    % Remove all f0 > fCut

[F00,ix]      = sort(Fxk(1,:));
nxk           = sum(~isnan(F00));
   
ixBest = 0;
for j=1:nxk
  i           = ix(j);
  NLP_x       = x_k(:,i);
  NLP_f       = Fxk(1,i);
  ProbL.x_0   = NLP_x;
  ProbL.f_0   = NLP_f;
  if NLP_f >= fCut
     Fxk(2,i) = inf;
  else
    if ~isempty(IntVars)
       x00 = max(x_L(IntVars),min(x_U(IntVars),ProbL.x_0(IntVars)));
       ProbL.x_L(IntVars)=x00;
       ProbL.x_U(IntVars)=x00;
    end
    rl       = tomRunFast(localSolver,ProbL);
    PrintResult(rl,tomPrint);
    totfEval = totfEval + rl.FuncEv;
    totcEval = totcEval + rl.ConstrEv;
    gradEv   = gradEv + rl.GradEv;
    HessEv   = HessEv + rl.HessEv;
    if M > 0
       % Check if point feasible enough to be accepted
       h_L1  = consviolation(rl, 0);
       if h_L1 <= M*ProbL.optParam.cTol
          f_k  = rl.f_k;
          if PriLev > 2
             fprintf('   ')
             fprintf('Local point %3d accepted, f(x) = %20.10f.',i,f_k);
             fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
             fprintf(' sum(|constr|)%18.10f',h_L1);
             if isempty(IntVars)
                fprintf('\n');
             else
                xprinti(rl.x_k(IntVars),' IV:');
             end
          end
       else
          f_k = Inf;
          if PriLev > 3
             fprintf('   ')
             fprintf('Local point %3d NOT OK,   f(x) = %20.10f.',i,rl.f_k);
             fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
             fprintf(' sum(|constr|)%18.10f',h_L1);
             if isempty(IntVars)
                fprintf('\n');
             else
                xprinti(rl.x_k(IntVars),' IV:');
             end
          end
       end
    else
       f_k = rl.f_k;
       if PriLev > 2
          fprintf('   ')
          fprintf('Local point %3d, f(x) = %20.10f.',i,f_k);
          fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
          if isempty(IntVars)
             fprintf('\n');
          else
             xprinti(rl.x_k(IntVars),' IV:');
          end
       end
    end
    Fxk(2,i) = f_k;
    XK(:,i)  = rl.x_k;
    if (fMinEQ > 0 & ~isinf(f_k) & isinf(fMin)) | f_k <= fMin
       fMin0     = fMin;
       if f_k == fMin
          xMin   = [xMin,rl.x_k];
          cDev   = [cDev,rl.c_k];
          ixBest = [ixBest,i];
       else
          xMin   = rl.x_k;
          cDev   = rl.c_k;
          ixBest = i;
       end
       fMin      = f_k;
       fMinEQ    = 0; % set better fMinEQ later on
       CMin      = ProbL.x_L;
       CMin(ixD) = (rl.x_k(ixD) - x_L(ixD))./x_D(ixD);
       if ~isempty(IntVars)
          IMin   = ProbL.x_L(IntVars);
       end
       if PriLev > 1 
          fprintf('   ')
          fprintf('Step %3d. Pnt %3d. ',j,i)
          if fMin0 == glcfMin
             fprintf('Best from Phase 2  ')
          else
             fprintf('Previous Best ')
          end
          if fMin0 >= 1E300
             fprintf('%18.16f.     ',inf)
          elseif fMin0 > 1E15
             fprintf('%18.16e.',fMin0)
          else
             fprintf('%18.16f.     ',fMin0)
          end
          fprintf(' New Best %18.16f',f_k)
          if isempty(IntVars)
             fprintf('\n');
          else
             xprinti(rl.x_k(IntVars),'. IV:');
          end
       end
    end
    if ~isempty(IntVars)
       ProbL.x_L(IntVars)=x_L(IntVars);
       ProbL.x_U(IntVars)=x_U(IntVars);
    end
    if PriLev > 5 
       fprintf('   ')
       fprintf('End of local search # %d out of %d\n',i,size(x_k,2)); 
       fprintf('                      ')
       fprintf('----------------------------\n\n'); 
    end
    if totfEval > mFunc, break; end
  end
end
if fMinEQ == 1000000, fMinEQ = 0; end % Reset if fMinEQ used as flag

if ixBest(1) > 0
   if PriLev > 0 
      fprintf('glcCluster Ph-3: ')
      fprintf('Local try')
      fprintf(' %d',ixBest)
      fprintf('/%d: ',nxk)
      fprintf('fMin %18.16f.',fMin)
      fprintf('  Ph-2 glcfMin ')
      if glcfMin >= 1E300
         fprintf('%18.16f.',Inf)
      else
         fprintf('%18.16f.',glcfMin)
      end
      if isempty(IntVars)
         fprintf('\n');
      else
         xprinti(xMin(IntVars),' IV:');
         for j=2:size(xMin,2)
             xprinti(xMin(IntVars,j),' IV:');
         end
      end
   end
else
   if PriLev > 0
      fprintf('glcCluster Ph-3: ')
      fprintf('No improvement in %d local tries.\n',nxk)
  end
end
% Update X0, XX and FX0
if nxk > 0
   X0  = [X0,x_k(:,ix(1:nxk))];
   XX  = [XX,XK(:,ix(1:nxk))];
   FX0 = [FX0,Fxk(:,ix(1:nxk))];
   nX0 = size(X0,2);
end

Cd=C;
if isempty(maxFunc2)
   if DFREE
      if mFunc < 3000
         maxFunc = max(n, mFunc - totfEval - min(200,max(80,20*n)));
      else
         %maxFunc = max(n, mFunc - totfEval - min(1000,max(400,100*n)));
         maxFunc = max(n, mFunc - totfEval - min(2000,max(800,200*n)));
         %maxFunc = max(n, mFunc - totfEval - 3000);
      end
   else
      %HKH maybe change this
      maxFunc = min(maxFunc1,max(n, mFunc - totfEval - 5*n));
   end
else
   maxFunc = maxFunc2;
end
Prob.optParam.MaxFunc = maxFunc;
Prob.optParam.MaxIter = maxFunc;
%Prob.optParam.MaxFunc = 2000;
%Prob.optParam.MaxIter = 2000;

if length(ignoreIdx) == length(F) | ~isempty(IntVars)
   % Skip Phase 4/5/6 if Integer problem with several choices
   Phase4 = 0;
else
   % ==========================================
   % Phase 4 Second DIRECT run
   % ==========================================
   % tomPrint                =2
   % Prob.optParam.IterPrint = 1;
   % Prob.PriLevOpt          = 2;
   Phase4                    = 1;
   Prob.WarmStart            = 1;
   Prob.optParam.IterPrint   = PriLev > 5;
   % Set best point in Prob for input to DIRECT solver
   if ~isinf(fMin)
      Prob.MIP.fIP              = fMin;
      Prob.MIP.xIP              = xMin(:,end);
   end
   glcfMin0                  = glcfMin;  % Save best DIRECT f(x) for test later
   fMinEQ0                   = fMinEQ;   
   Result                    = tomRunFast(DIRECT,Prob);
   PrintResult(Result,tomPrint);
   totfEval                  = totfEval + Result.FuncEv;
   totcEval                  = totcEval + Result.ConstrEv;

   load(DIRECTFILE,'feasible'); 
   REJECT = 0;
   if M > 0 & feasible
      if isempty(Result.x_k)
         % No solution found
      else
          Result.c_k = nlp_c(Result.x_k(:,1),Prob);
          % Check if point feasible enough to be accepted
          h_L1 = consviolation(Result, PriLev-1);
          if h_L1 > 10*M*ProbL.optParam.cTol
             % Infeasible second glcDirect/glcSolve solution
             REJECT = 1;
          end
       end
    end

   % ==========================================
   % Phase 5 Second clustering
   % ==========================================

   % Compute clustering for full set of DIRECT points
   % Should be no integer variables
   if DFREE & (mFunc-totfEval) < 300
      minLocalTry               = 3; 
   else
      minLocalTry               = 6; 
   end
   maxLocalTry               = 30; 

   % C and F could have been returned directly from glcDirect/glcSolve
   if OLD
      load(DIRECTFILE,'C','glcfMin','F','ignoreIdx','fMinIdx','fMinEQ','G'); 
   else
      load(DIRECTFILE,'C','glcfMin','F','ign','fMinIdx','fMinEQ','G'); 
      ignoreIdx = ign;
   end
   if PriLev > 0 
      fprintf('glcCluster Ph-4: ')
      fprintf('Tried %4d. ',Result.FuncEv)
      if REJECT
         fprintf('New global best f ')
         fprintf('%16.10f (L1 infeasibility = %f)\n',Result.f_k,h_L1)
      elseif glcfMin < fMin & fMinEQ == 0
         fprintf('New global best f ')
         fprintf('%20.10f\n',glcfMin)
      elseif glcfMin0 <= glcfMin
         % No improvement in Phase 4 compared to Phase 1 by glcDirect/glcSolve
         fprintf('No improvement. ')
         if fMinEQ > 0
            fprintf('Old fMinEQ %12.10f. New fMinEQ %12.10f',fMinEQ0,fMinEQ);
         end
         fprintf('\n')
      else
         fprintf('Improved global solution glcfMin ')
         fprintf('%20.10f. ',glcfMin)
         if fMinEQ > 0
            fprintf('Old fMinEQ %12.10f. New fMinEQ %12.10f',fMinEQ0,fMinEQ);
         end
         fprintf('\n')
      end
   end
   % Update fMin if glcfMin from DIRECT in Phase 4 better than Phase 3 solution
   if glcfMin < fMin & ~REJECT & fMinEQ == 0
      fMin = glcfMin;
      xMin = Result.x_k;
      cDev = Result.c_k;
   end
   if M == 0
      iPnt  = ones(1,length(F));
      ix    = [1:length(F)];
   else
      iPnt  = double(all(G <= cTol));
      ix    = find(iPnt > 0);
      if isempty(ix)
         cTol = 100*cTol;
         if PriLev > 0
            fprintf('   Increasing cTol %f\n',cTol);
         end
         iPnt  = double(all(G <= cTol));
         ix    = find(iPnt > 0);
      end
      if isempty(ix)
         % Desperate - avoid crash
         iPnt  = ones(1,length(F));
         ix    = [1:length(F)];
      end
      % Must set nonused elements to <0, because integer keys from 0,1,...
      iPnt(iPnt==0) = -inf;
   end


   F0        = [];
   n0        = 0;
   c_0       = [];
   nOK       = 0;
   minlV0    = [];
   v2V0      = [];
   maxDistV0 = [];

   iVal     = 1;
   nVal     = 1;
   Fout     = Inf*ones(nVal,1);
   Cout     = zeros(n,nVal);

   % Update points to ignore, not considered for trisection 
   % Later all points that are clustered are set to be ignored
   ignore            = zeros(length(F),1);
   ignore(ignoreIdx) = 1;
   % Best single point sampled, not included in any cluster
   nPnts             = min(1000,maxLocalTry*nVal);
   Fv                = zeros(nPnts,1);
   c_k               = zeros(n,nPnts);
   nv                = 0;
   maxDistV          = NaN*ones(nVal,1);
   maxDist           = NaN;
   minlV             = [];
   v2V               = [];
   i0                = nOK;
   for i = 1:nVal
      % Indices with fixed integer combination corresponding to key iVal(i)
      iKey = iVal(i);
      ix   = find(iPnt == iKey);
      % Call cluster algorithm
      [FvNew,cNew,Fout(i0+i),Cout(Reals,i0+i),maxDist,minl,v2] = clustAlg(...
         C(Reals,ix),F(ix),PriLevCl,minLocalTry,maxLocalTry,incLocalTry,...
         incDist,decDist,maxDistMin);
   
      nPnt                    = length(FvNew);
      Fv(nv+1:nv+nPnt)        = FvNew;
      if ~isempty(IntVars)
         c_k(IntVars,nv+1:nv+nPnt) = C(IntVars,ix(1:nPnt));
      end
      c_k(Reals,nv+1:nv+nPnt) = cNew;
      nv                      = nv + nPnt;
      maxDistV(i)             = maxDist;
      if isempty(minlV)
         minlV         = minl;
      else
         minlV         = [minlV,0,minl];
      end
      if isempty(v2V)
         v2V         = v2;
      else
         v2V         = [v2V,0,v2];
      end
   end
   
   % Safer to always use the best far away point
   ix = [1:nOK+nVal];
   nx = length(ix);
   if nx > 0
      % Add start from best non-clustered points
      Fv(nv+1:nv+nx)    = Fout(ix);
      c_k(:,nv+1:nv+nx) = Cout(:,ix);
   end
   nPnt = nv + nx;
   if n0 > 0
      % Add start from other infeasible integer combinations
      Fv(nPnt+1:nPnt+n0)    = F0;
      c_k(:,nPnt+1:nPnt+n0) = c_0;
      nPnt                  = nPnt + n0;
   end
   if nPnt < nPnts
      % Shrink to used size
      Fv  = Fv(1:nPnt);
      c_k = c_k(:,1:nPnt);
   end
   if PriLev > 4
      xprint(Fv,'Fv: ')
      xprint(c_k,'c_k:',[],n)
   end
   nTot         = nPnt;
   [c_k I J]    = unique(c_k','rows');
   c_k          = c_k';
   Fv           = Fv(I); 
   nPnt         = size(c_k,2);
   x_k          = tomsol(9,x_L,c_k,x_D);
   x_k(IntVars) = round(x_k(IntVars));

   if PriLev > 0
      fprintf('   ')
      fprintf('Total pnts  %d. ',nTot)
      fprintf('Cluster pnts %d. ',nv)
      fprintf('Single  pnts %d. ',nVal)
      if ~isempty(IntVars)
         fprintf('Infeasible cluster pnts %d.',n0)
         fprintf('Single nonfeasible pnts %d. ',nOK)
      end
      if nPnt ~= nTot
         fprintf(' Unique pnts %d!',nPnt)
      end
      fprintf('\n')
   end

   nxk = size(x_k,2);
   if PriLev > 0
      fprintf('glcCluster Ph-5: ')
      fprintf('%d initial points from clustering.',nxk)
      fprintf('\n')
   end
   Fxk = zeros(2,nxk);
   XK  = zeros(n,nxk);
   if 0 & fMinEQ > 0
      for i = 1:nxk
          Fxk(1,i) = nlp_f(x_k(:,i), ProbL);         % Compute f(x_0) for each x_0
          totfEval = totfEval + 1;
      end 
   else
      Fxk(1,:) = Fv';
   end

   % Remove same/similar initial values
   % HKH New Maybe check FX0(2,:), i.e. no start in old optimum
   remove = 0;
   for P = 1:nX0
       f0   = FX0(1,P);
       if ~isnan(f0)
          ix   = find(abs(f0-Fxk(1,:)) < fEqTol*max(1,abs(f0)));
          if ~isempty(ix)
             iE              = checkEq(X0(:,P),x_k(:,ix),xEqTol);
             if ~isempty(iE)
                remove       = remove + length(iE);
                % HKH New
                %F0(1,ix(iE)) = NaN;
                Fxk(1,ix(iE)) = NaN;
             end
          end
       end
   end
   Fxk(1,Fxk(1,:) >= fCut) = NaN;                    % Remove all f0 > fCut

   [F00,ix]      = sort(Fxk(1,:));
   nxk           = sum(~isnan(F00));

end

% *******************************************************************
if ~Phase4
   Phase6 = 0;
elseif isempty(maxFunc2)
   Phase6 = (mFunc - totfEval) > 0;
else
   Phase6 = maxFunc2 > 0;
end

% ==============================================
% Phase 6 Local search from best cluster points
% ==============================================
if Phase6
   ProbL.optParam.MaxFunc  = max(5*ProbL.N,mFunc-totfEval);
   %ProbL.optParam.MaxIter = ProbL.optParam.MaxFunc;
   ProbL.optParam.MaxIter  = mIter;
   ProbL.WarmStart         = 0;
   ixBest                  = [];
   for j=1:min(100,nxk)
       i           = ix(j);
       NLP_x       = x_k(:,i);
       NLP_f       = Fxk(1,i);
       ProbL.x_0   = NLP_x;
       ProbL.f_0   = NLP_f;
      if NLP_f >= fCut
         Fxk(2,i) = inf;
      else
       if ~isempty(IntVars)
          x00 = max(x_L(IntVars),min(x_U(IntVars),ProbL.x_0(IntVars)));
          ProbL.x_L(IntVars)=x00;
          ProbL.x_U(IntVars)=x00;
       end
       rl   = tomRunFast(localSolver,ProbL);
       PrintResult(rl,tomPrint);
       totfEval=totfEval + rl.FuncEv;
       totcEval=totcEval + rl.ConstrEv;
       gradEv=gradEv + rl.GradEv;
       HessEv=HessEv + rl.HessEv;
       if M > 0
          % Check if point feasible enough to be accepted
          h_L1 = consviolation(rl, PriLev-1);
          if h_L1 <= M*ProbL.optParam.cTol
             f_k = rl.f_k;
             if PriLev > 1
                fprintf('   ')
                fprintf('Local point accepted, f(x) = %20.10f.',f_k);
                fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
                if isempty(IntVars)
                   fprintf('\n');
                else
                   xprinti(rl.x_k(IntVars),' IV:');
                end
             end
          else
             f_k = Inf;
          end
       else
          f_k = rl.f_k;
          if PriLev > 2
             fprintf('   ')
             fprintf('Local point %3d, f(x) = %20.10f.',i,f_k);
             fprintf(' f(x_0) = %20.10f.',ProbL.f_0);
             if isempty(IntVars)
                fprintf('\n');
             else
                xprinti(rl.x_k(IntVars),' IV:');
             end
          end
       end
       Fxk(2,i) = f_k;
       XK(:,i)  = rl.x_k;
       if f_k <= fMin
          % Local point does improve
          fMin0 = fMin;
          if f_k == fMin
             xMin   = [xMin,rl.x_k];
             cDev   = [cDev,rl.c_k];
             ixBest = [ixBest,i];
          else
             xMin   = rl.x_k;
             cDev   = rl.c_k;
             ixBest = i;
          end
          fMin = f_k;
          CMin = ProbL.x_L;
          CMin(ixD) = (xMin(ixD) - x_L(ixD))./x_D(ixD);
          if ~isempty(IntVars)
             IMin = ProbL.x_L(IntVars);
          end
          if PriLev > 1 
             fprintf('   ')
             fprintf('Step %d. Pnt %d. ',j,i)
             fprintf('Best Ph-4 %18.16f. ',glcfMin)
             fprintf('Old Best %18.16f. ',fMin0)
             fprintf('New Best %18.16f',f_k)
             fprintf('\n')
          end
       end
       if ~isempty(IntVars)
          ProbL.x_L(IntVars)=x_L(IntVars);
          ProbL.x_U(IntVars)=x_U(IntVars);
       end
       if totfEval > mFunc, break; end
      end
   end
   % Update X0, XX and FX0
   if nxk > 0  % Use j <= nxk
      X0  = [X0,x_k(:,ix(1:j))];
      XX  = [XX,XK(:,ix(1:j))];
      FX0 = [FX0,Fxk(:,ix(1:j))];
      nX0 = size(X0,2);
   end
end
if Phase6 & fMin < glcfMin
   if PriLev > 0 
      fprintf('glcCluster Ph-5: ')
      if isempty(ixBest)
         fprintf('Tried %d out of %d. No improvement. ',j,nxk)
         fprintf('fMin: %20.10f ',fMin)
      else
         fprintf('Best try')
         fprintf(' %d',ixBest)
         fprintf('. Tried %d out of %d. ',j,nxk)
         fprintf('fMin: %20.10f ',fMin)
      end
      if glcfMin >= 1E300
         fprintf('glcfMin: %20.10f\n',Inf)
      else
         fprintf('glcfMin: %20.10f\n',glcfMin)
      end
   end
   feasible   = 1;
   % Set best point in Prob for input to DIRECT solvers
   if ~isinf(fMin)
      Prob.MIP.fIP              = fMin;
      Prob.MIP.xIP              = xMin(:,end);
   end
elseif Phase6
   if PriLev > 0 
      fprintf('glcCluster Ph-5: ')
      fprintf('No improvement in %d local tries. ',i)
      fprintf('fMin: %20.10f ',fMin)
      if glcfMin >= 1E300
         fprintf('Best Global: %20.10f\n',Inf)
      else
         fprintf('Best Global: %20.10f\n',glcfMin)
      end
   end
end

if fMin >= 1E300 & fMinIdx > 0
   % Safe guard if no feasible point found
   load(DIRECTFILE,'F');
   fMin = F(fMinIdx);
end

% ********************************************************************

Result.GradEv   = gradEv;
Result.HessEv   = HessEv;
Result.x_k      = xMin;        % x corresponding to best f
Result.f_k      = fMin;        % Best function value
Result.c_k      = cDev;        % Constraint values at best function value
Result.Ax       = [];          % Reset as empty for safety
Result.Iter     = 1;           % Number of iterations
Result.FuncEv   = totfEval;
Result.ConstrEv = totcEval;
if feasible
   Result.ExitFlag = 0;
else
   Result.ExitFlag = 7;
end
Result.ExitText = ['Tried ' num2str(totfEval+fPrev) ...
                   ' function values in total '];
Result.Solver          = 'glcCluster';
Result.SolverAlgorithm = ...
 ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];

% Save cluster information
Result.Cluster.maxDist = maxDistV;
Result.Cluster.minDist = minlV;
Result.Cluster.minIdx  = v2V;
Result.Cluster.x_k     = x_k;
Result.Cluster.f_k     = Fv;
Result.Cluster.c_k     = c_k;

% Save local tries
Result.Cluster.FX0     = FX0;
Result.Cluster.X0      = X0;
Result.Cluster.XX      = XX;
Prob.GO.ProbL          = ProbL;
Result                 = endSolve(Prob,Result);

if PriLev > 10 
   fprintf('***********************************************************\n');
   totfEval
   gradEv
   Result.x_k
   fprintf('***********************************************************\n\n\n');
end

%----------------------------------------------------------------------------
function List=makeAllNorms(X)
%----------------------------------------------------------------------------
n = size(X,2);

%List1 = 100*ones(1,n);

%for i=1:size(X,2)
%    minX = 100;
%    for j=1:i-1
%        minX = min(minX,norm(X(:,i)-X(:,j)));
%    end
%    for j=i+1:n
%        minX = min(minX,norm(X(:,i)-X(:,j)));
%    end
%    List1(i) = minX;
%end

%%D = zeros(n,n);
%D = sparse([],[],[],n,n,n*(n+1)/2);
%for i=1:n-1
%    D(i,i+1:n) = sum((X(:,i)*ones(1,n-i)-X(:,i+1:n)).^2);
%end
%
%List = sqrt(min(D+D'+1000*speye(n,n)));

%D=D+D'+1E10*eye(n,n);
%List = sqrt(min(D));

% Make more memory efficient, avoid square root

z = sum((X(:,1)*ones(1,n-1)-X(:,2:n)).^2);
List = [min(z),z];
for i=1:n-1
    z = sum((X(:,i)*ones(1,n-i)-X(:,i+1:n)).^2);
    List(i) = min([List(i),z]);
    List(i+1:n) = min(List(i+1:n),z);
end

%----------------------------------------------------------------------------
function clusters=cluster(X,maxDist)
%----------------------------------------------------------------------------

clusters=[];

[md idx]=min(sqrt(sum(X.^2,1)));
%[md idx]=min(sum(X.^2,1));

idx=idx(1); % We are only interested in one of these points.
  % If the other points are in the same cluster they will be found any way.
  % If they should start there own clusters, they will be found later on.

            
clusters(1) = idx;         
clCnt = 1;
ok = 1;
k  = idx;

while ok
   tmpIdx=ones(1,size(X,2));
   tmpIdx(clusters(clusters~=0)) = 0;
   S = find(tmpIdx==1);
   sLen = length(S);
   if sLen == 0
      ok = 0;
   else
      %distVec = sqrt(sum((X(:,k)*ones(1,sLen) - X(:,S)).^2,1));

      distVec=tomsol(30,X(:,k),X(:,S));

      idx = find(distVec<=maxDist);

      if ~isempty(idx)
         clusters = [clusters S(idx)];
      else
         if clCnt == length(clusters)
            [md idx]= min(distVec);
            clusters = [clusters 0 S(idx(1))];
            clCnt = length(clusters)-1;
         end
      end

      clCnt = clCnt + 1;
      k = clusters(clCnt);
      if k == 0
         clCnt = clCnt + 1;
         k = clusters(clCnt);
      end
   end
end
%----------------------------------------------------------------------------
function [Fv,c_k,Fout,Cout,maxDist,minl,v2] = clustAlg(C,F,PriLev,minClPnt,...
                           maxClPnt,incClPnt,incDist,decDist,maxDistMin)
%----------------------------------------------------------------------------

% makeAllNorms returns the distance from any point
% to that points closest neighbour.
% minL2(i) = min (| C(:,i) - C(:,j)|, j=1:m, j~=i, m=size(C,2))
% minL2    = makeAllNorms(C);

if size(C,2) == 1
   meanDist = inf;
   stopDist = inf;
   v1       = 1;
   v2       = [];
else
   minL2    = tomsol(29,C);

   % meanDist=mean( minL2 ) is the average distance any point has to its closest neighbour
   meanDist = mean(minL2);
   stopDist = meanDist + 0.9*(max(minL2)-meanDist);
   % meanDist will be an upper bound on maxDist
   % Points with distance >= meanDist wont be clustered, vIdx == 1
   vIdx     = minL2 >= (meanDist); % Safeguard if all values are very close
   
   % v2 is the set of points being considered for clustering, possibly []
   v2       = find(vIdx == 0);

   %v2Len = length(v2)
end

if isempty(v2) | length(v2) == size(C,2)
   % All points are on the same distance to each other
   % Use lowest F in (Cout,Fout) and 
   % up to maxClPnt points with lowest F in (Fv,c_k)
   if length(v2) == size(C,2)
      v1 = v2;
      v2 = [];
   end

   [FF v1] = sort(F(:)');
   Fout    = F(v1(1));
   Cout    = C(:,v1(1));
   nPnt    = min(maxClPnt,length(F)-1);
   v2      = v1(2:nPnt+1); 
   Fv      = F(v2);
   c_k     = C(:,v2);
   maxDist = meanDist;
   minl    = v1(1:nPnt);
   if PriLev > 1
      fprintf('     1:');
      fprintf(' clustAlg: ');
      fprintf('EquDist %8.5f. ',maxDist)
      fprintf('# of clusters %3d. ',nPnt)
      fprintf('Points used %5d. ',length(v2))
      fprintf('Total points %d. ',size(C,2))
      fprintf('\n')
   end
else
   % Best point among the non clustered points
   v1          = find(vIdx);
   %v1Len = length(v1)
   [Fout Iout] = min(F(v1));
   Cout        = C(:,v1(Iout));
   % Points that are to be clustered will no longer be considered for trisection
   ignore(v2)        = 1;

   % Now, for this subset, find the closest point to all point
   minl=tomsol(29,C(:,v2));

   % Perhaps this min number should be changed so that 
   %         if max(minl) < 0.01 maxDist = 0.05
   %         if  0.01 < max(minl) < 0.1 maxDist = 0.1
   %         else maxDist = max(minl)
   if 1
      maxDist = max([max(minl) maxDistMin]);
      % Avoid square root
      %maxDist = max([max(minl) 0.05^2]);
   else
      maxminl=max(minl);
      if maxminl < 0.01 
         maxDist = 0.05;
      elseif 0.01 < maxminl < 0.1 
         maxDist = 0.1;
      else 
         maxDist = maxminl;
      end
   end

   v3   = cluster(C(:,v2),maxDist);
   v4   = find(v3==0);
   nPnt = 1+length(v4);

   if PriLev > 1
      fprintf('     2:');
      fprintf(' clustAlg: ');
      fprintf('maxDist %8.5f. ',maxDist)
      fprintf('# of clusters %3d. ',nPnt)
      fprintf('Points used %5d. ',length(v2))
      fprintf('Total points %d. ',size(C,2))
      fprintf('\n')
   end

   if nPnt > maxClPnt
      %while nPnt > maxClPnt & maxDist < meanDist
      while nPnt > maxClPnt & maxDist < stopDist
         % Increase maxDist
         maxDist = incDist*maxDist;
         v3  = cluster(C(:,v2),maxDist);
         v4  = find(v3==0);
         nPnt = 1+length(v4);
         if PriLev > 1
            fprintf('     3:');
            fprintf(' clustAlg: ');
            fprintf('maxDist %8.5f. ',maxDist)
            fprintf('# of clusters %3d. ',nPnt)
            fprintf('Points used %5d. ',length(v2))
            fprintf('Total points %d. ',size(C,2))
            fprintf('\n')
         end
      end
   elseif nPnt < minClPnt
      maxDist = decDist*maxDist;
      while nPnt < minClPnt & maxDist > maxDistMin
         v3S = v3;
         v4S = v4;
         v3  = cluster(C(:,v2),maxDist);
         v4  = find(v3==0);
         nPnt = 1+length(v4);
         if PriLev > 1
            fprintf('     4:');
            fprintf(' clustAlg: ');
            fprintf('maxDist %8.5f. ',maxDist)
            fprintf('# of clusters %3d. ',nPnt)
            fprintf('Points used %5d. ',length(v2))
            fprintf('Total points %d. ',size(C,2))
            fprintf('\n')
         end
         maxDist = decDist*maxDist;
      end
      if nPnt > maxClPnt
         % Too many clusters, use previous clustering
         v3      = v3S;
         v4      = v4S;
         maxDist = maxDist/decDist;
         nPnt    = 1+length(v4);
         if PriLev > 1
            fprintf('     5:');
            fprintf(' clustAlg: ');
            fprintf('Prev. clustering, ')
            fprintf('# of clusters %3d. ',nPnt)
            fprintf('Points used %5d. ',length(v2))
            fprintf('Total points %d. ',size(C,2))
            fprintf('\n')
         end
      end
      maxDist = maxDist/decDist;
   elseif nPnt < maxClPnt & incClPnt
      maxDist = decDist*maxDist;
      while nPnt < maxClPnt & maxDist > maxDistMin
         v3S = v3;
         v4S = v4;
         v3  = cluster(C(:,v2),maxDist);
         v4  = find(v3==0);
         nPnt = 1+length(v4);
         if PriLev > 1
            fprintf('     6:');
            fprintf(' clustAlg: ');
            fprintf('maxDist %8.5f. ',maxDist)
            fprintf('# of clusters %3d. ',nPnt)
            fprintf('Points used %5d. ',length(v2))
            fprintf('Total points %d. ',size(C,2))
            fprintf('\n')
         end
         maxDist = decDist*maxDist;
      end
      if nPnt > maxClPnt
         % Too many clusters, use previous clustering
         v3      = v3S;
         v4      = v4S;
         maxDist = maxDist/decDist;
         nPnt    = 1+length(v4);
         if PriLev > 1
            fprintf('     7:');
            fprintf(' clustAlg: ');
            fprintf('Prev. clustering, ')
            fprintf('# of clusters %3d. ',nPnt)
            fprintf('Points used %5d. ',length(v2))
            fprintf('Total points %d. ',size(C,2))
            fprintf('\n')
         end
      end
      maxDist = maxDist/decDist;
   end
   nv4 = length(v4);
   if nv4 >= 1
      ix            = v2(v3(1:v4(1)-1));
      [Fv(1) idxMt] = min(F(ix));
      c_k(:,1)      = C(:,ix(idxMt(1)));
      for i=2:nv4
          ix            = v2(v3(v4(i-1)+1:v4(i)-1) );
          [Fv(i) idxMt] = min(F(ix));
          c_k(:,i)      = C(:,ix(idxMt(1)));
      end
      ix                = v2(v3(v4(nv4)+1:length(v3)));
      [Fv(1+nv4) idxMt] = min(F(ix));
      c_k(:,1+nv4)      = C(:,ix(idxMt(1)));
   else
      ix            = v2(v3);
      [Fv idxMt]    = min(F(ix));
      c_k(:,1)      = C(:,ix(idxMt(1)));
   end
end
if PriLev >= 1
   fprintf('   CLUSTALG: ');
   fprintf('     maxDist %8.5f. ',maxDist)
   fprintf('# of clusters %3d. ',length(Fv))
   fprintf('Points used %5d. ',length(v2))
   fprintf('Total points %d. ',size(C,2))
   fprintf('meanDist %8.5f. ',meanDist)
   fprintf('stopDist %8.5f. ',stopDist)
   fprintf('\n')
end
if isempty(Fout), keyboard, end
if isempty(Cout), keyboard, end

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

   
% MODIFICATION LOG:
%
% 011103  hkh  Written from see test algorithm
% 011106  hkh  Rewriting cluster and makeAllNorms more efficiently
% 011110  hkh  Fixed errors in comments
% 011112  hkh  Improvement in the use of NPSOL and SNOPT (internal derivatives)
% 011113  hkh  Use Prob.GO for maxFunc1 and maxFunc2, besides additional input
%              Use any derivative information, use internal SOL derivatives
% 020110  hkh  Change name fMin to glcfMin. Conflict in 5.x with fmin function
% 020110  hkh  Loop to get feasible in 1st step.
% 020111  hkh  Safe guard maxFunc1 to odd number. Print improved function values
% 020411  hkh  Speedup bottle necks (distance computations) using Fortran code
% 020412  hkh  Setting fMin, not glcfMin and fMinIdx missing
% 020412  hkh  Avoid division with 0, if fixed variables
% 020701  hkh  Call GetSolver for default local solver
% 031126  hkh  Handle integer variables by fixing them in local search
% 031201  hkh  Avoid doing preSolve. Changed glcSolve to glcCluster in output
% 031204  hkh  Use additional structure ProbL for local search
% 031204  hkh  Check constraint violation of local solution before accepting it
% 031204  hkh  Wrong index when putting back best local solution to mat file
% 031204  hkh  Major revision, bug fixes
% 040110  hkh  xMin not set, if local searched failed to improve
% 040111  hkh  Accept local point, even if glcFast: f_k + fMinEQ is lower
% 040203  hkh  Expression idx3 = find(all( C==c_k(:,i)*ones(1,size(C,2)) )); 
%              does not work for n==1. If n==1 remove all( )
% 040308  hkh  Give ERROR if called with pure IP
% 040308  hkh  New input GO.DIRECT, change glcFast to glcSolve optional
% 040308  hkh  Must do ixMin = fMinIdx(1); if many solutions
% 040326  hkh  Must test Phase IV result even if no Phase6
% 040329  hkh  Wrong initial points selected in Phase V
% 040329  hkh  Save equal minima in Phase V, and return > 1 global minima
% 040402  hkh  New option WarmStart, full call to iniSolve
% 040404  hkh  Major revision
% 040604  hkh  Add output field Cluster, and input localTry and maxDistMin
% 040604  med  Added sorting of cluster points, x_k and Fv
% 040609  hkh  Improve comments, more printout before keyboard command
% 040928  hkh  Print localSolver in SolverAlgorithm
% 041006  hkh  Revision of integer handling
% 041017  hkh  Use MIP.fIP and MIP.xIP to send best point to glcFast/Solve
% 041018  hkh  Print number of trial points in Phase 4
% 041018  hkh  cDev only updated for Phase 4 if point accepted as feasible
% 041123  hkh  Change call to tomRun
% 041129  hkh  Change field GO.Prob to GO.ProbL in both Prob and Result
% 041129  hkh  Prob.GO.ProbL=ProbL, implies Result.Prob.GO.ProbL always output 
% 041129  hkh  2nd DIRECT call, default < maxFunc1 function evals
% 041129  hkh  Change comments for maxFunc2 (restriced by maxFunc1 by default)
% 041129  hkh  Change comments for maxDist, minDist, maxDistMin
% 050220  hkh  maxFunc3 as input was not OK treated
% 050220  hkh  maxFunc3 is Prob.GO.maxFunc3, or input maxFunc3, or 3*maxFunc1
% 060327  hkh  Major revision of comments and avoid unnecessary copying of arrays
% 060327  hkh  Add name of DIRECT solver in input
% 060327  hkh  Added counting on number of constraint evaluations, Result.ConstrEv
% 060327  hkh  Change default to glcDirect, glcFast still possible to use
% 060327  hkh  Set Result.Ax=[], otherwise wrong value might be returned
% 060327  hkh  Always test if to utilize known derivatives, use internal SOL derivs
% 060327  hkh  Always reset the optParam structure for the local solver
% 060327  hkh  Implement Adaptive Clustering procedure, new GO inputs
% 060328  hkh  Skip keyboard, and display error if bad maxFunc input
% 060814  med  FUNCS used for callbacks instead
% 070222  hkh  Revise IntVars handling, use new format
% 070222  hkh  Change output in step loop, print 3 different f(x) values
% 070303  hkh  Greatly improved algorithm for integer and constrained problems
% 070305  hkh  Revised code. New sub function clustAlg, prepare for separate code
% 070307  hkh  clustAlg: if either v1 or v2 [], sort(F) and take best points
% 070531  hkh  clustAlg: sort(F) gave col result, crash, change to sort(F(:)')
% 070921  hkh  ix not defined if M==0
% 071009  hkh  Add local searches from all points in Prob.x_0
% 071010  hkh  Avoid Prob.x_0=0 initial point set by Tomlab, if Prob.X0 defined
% 080205  hkh  nPnts set too big
% 080212  hkh  Always nPnts>1. Use stopDist >> meanDist to decrease clusterpnts
% 080229  hkh  Set Prob.solvType after ProbL defined
% 080416  hkh  Avoid taking advantages of derivatives - could be dummy routines
% 080417  hkh  DIRECT used relative violation, consviolation absolute==> bad solution
% 080417  hkh  Now relative consviolation makes possible small local improvements
% 091206  hkh  Major revision, major algorithm changes
% 100217  hkh  Error twice comparing (X0,FX0) with (x_k,Fxk), avoiding same x_0
% 100714  hkh  Avoid creating file SOL.txt with output for DFREE == 1
% 110615  hkh  Use fMin = Prob.MIP.fIP if WarmStart. Def field in WarmStartGLOBAL
% 110728  hkh  Matrix c_k was overwritten in Phase 6 loop
% 110728  hkh  Fixed warm start, adding fields from Result.Cluster into Prob.Cluster
