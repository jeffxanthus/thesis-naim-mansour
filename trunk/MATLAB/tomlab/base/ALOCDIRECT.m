% ALOCDIRECT.m
%
% ALOCDIRECT solves general constrained mixed integer global optimization
%
% ALOCDIRECT is a hybrid algorithm, that is using one of the following DIRECT
% algorithms: glcDirect (default),glcFast or glcSolve, for global search (Step 1). 
% Step 2 is an adaptive clustering algorithm to find a suitable number of clusters, 
% where the best point in each cluster is then used as an initial point for a 
% local search (Step 3).
% The 4th step is to run the DIRECT algoirithm once again, to possible improve.
% If the DIRECT algorithm improves the best point, a local search is 
% finally made as Step 5 with the new best point(s) as starting points.
% For a more detailed algorithm description, see below after USAGE:
%
% ALOCDIRECT solves problems of the form:
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
% function Result = ALOCDIRECT(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)
%
% INPUT PARAMETERS
%
% The following three parameters are either set in Prob.GO, or as extra input
%
% maxFunc1  Number of function evaluations in 1st call to the DIRECT solver.
%           Should be odd number (automatically corrected).
%           Default 100*dim(x) + 1.  May also be set as Prob.GO.maxFunc1
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
%           they are used. If only one structure Prob is given, ALOCDIRECT
%           uses the derivative routines given in the this structure
%           May also be set as Prob.GO.ProbL
%
% Prob    Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say GLCF
%   FUNCS.c   The routine to compute the nonlinear constraint, say GLCC
%             A call to tomFiles.m or glcAssign.m sets these fields.
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   PriLevOpt Print Level
%             0 = silent.
%             1 = Some output from each ALOCDIRECT phase
%             2 = More output from each phase
%             3 = Adaptive clustering information
%             6 = Use PrintResult( ,1) to print summary from each global and
%                 local run
%             7 = Use PrintResult( ,2) to print summary from each global and
%                 local run
%             8 = Use PrintResult( ,3) to print summary from each global and
%                 local run
%
%   WarmStart If true, >0, ALOCDIRECT warm starts the DIRECT solver.
%             The DIRECT solver will utilize all points sampled in last run,
%             from one or two calls, dependent on the success in last run.
%             Note: The DIRECT solver may not be changed if doing WarmStart
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
%             As 1st step needs maxFunc = 100*n + 1; maxFunc >> 100*n
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
% GO           Structure in Prob, Prob.GO.
% ---------------------------------------
%   localSolver Optionally change local solver used ('snopt' or 'npsol' etc.)
%   maxFunc1    See description above
%   maxFunc2    See description above
%   maxFunc3    See description above
%   ProbL       See description above
%
%   DIRECT      DIRECT subsolver, either glcDirect (default),glcFast or glcSolve
%   maxLocalTry Maximal number of local searches from cluster points
%               If <= 0, ALOCDIRECT stops after clustering. Default 30
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
%             Defines integer optimization parameters. Fields used:
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
%  ConstrEv Number of function evaluations
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
% Phase 1: Run DIRECT maxFunc1 function value trials
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
% Phase 3: Do a local search with each of the nPnt best cluster points as
%          initial starting value. Most likely the local search will find
%          the global optimum, if there are not too many local minima
% Phase 4: If the best point point in the local searches in Phase 3 is better 
%          than the best point found in the global search in Phase 1, this 
%          new point is added as input to the DIRECT solver and a warm
%          start run with DIRECT doing maxFunc2 function trials is done.
% Phase 5: If DIRECT has improved in Phase 4, a local search is done from
%          each of the points with the best function value
%          If the local search improves the best point found, then this point
%          is inserted in the DIRECT data base mat file used.
%          It is then possible for the user to do further warm start runs
%          with the DIRECT solver used 
%          A warm start of ALOCDIRECT is also possible.
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
%    Result = tomRun('ALOCDIRECT',Prob,2);
%
% Direct solver call:
%    Result = ALOCDIRECT(Prob);
% or with all input parameters as
%    Result = ALOCDIRECT(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)
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
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Nov 3, 2001.    Last modified Feb 22, 2007.

function Result = ALOCDIRECT(Prob, maxFunc1, maxFunc2, maxFunc3, ProbL)

if nargin < 5
   ProbL = [];
   if nargin < 4
      maxFunc3 = [];
      if nargin < 3
         maxFunc2 = [];
         if nargin < 2
            maxFunc1 = [];
            if nargin < 1
               error('ALOCDIRECT needs input structure Prob');
            end
         end
      end
   end
end

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
   disp('ALOCDIRECT requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'ALOCDIRECT requires both lower and upper variable bounds';
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
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
if IterPrint
   PriLev = max(1,PriLev);
end

if isinf(fGoal), fGoal = -1E300; end


if isfield(Prob.GO,'localSolver')
   localSolver   = deblank(Prob.GO.localSolver);
else
   localSolver = GetSolver('con',0,0);
end
if isfield(Prob.GO,'maxLocalTry')
   maxLocalTry   = Prob.GO.maxLocalTry;
else
   maxLocalTry = [];
end
if isempty(maxLocalTry), maxLocalTry = 30; end

if isfield(Prob.GO,'minLocalTry')
   minLocalTry   = Prob.GO.minLocalTry;
else
   minLocalTry = [];
end
if isempty(minLocalTry), minLocalTry = 1; end

Result                 = ResultDef(Prob);
Result.Solver          = 'ALOCDIRECT';
Result.SolverAlgorithm = ...
 ['Adaptive search with Constrained DIRECT solver ' DIRECT '-Local Search with ' localSolver];

% Do a preSolve, hopefully shrinking the box
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

Prob.x_L = x_L;
Prob.x_U = x_U;

n = Prob.N;

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
      error('ALOCDIRECT: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);
if length(IntVars) == n
   fprintf('ALOCDIRECT: Do not call ALOCDIRECT for PURE integer programs\n');
   fprintf('            ');
   fprintf('No point in clustering and local search for pure IP\n');
   fprintf('            ');
   fprintf('Call glcDirect, glcSolve, or some other IP solver instead\n');
   error('ALOCDIRECT: Wrong type of optimization problem');
end

if isempty(maxFunc1)
   if isfield(Prob.GO,'maxFunc1')
      maxFunc1 = Prob.GO.maxFunc1;
   end
end
if isempty(maxFunc1) 
   maxFunc1 = 100*n + 1;    % Odd number wanted :)
else
   % Make sure it is odd
   maxFunc1 = floor(maxFunc1/2)*2 + 1;
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
   %maxFunc = ceil(mFunc*0.8);
   %fprintf('Will now use %d for ',maxFunc);
   %fprintf('%s,',DIRECT); 
   %fprintf(' but a Warm Start might be needed!');
   %fprintf('\nNow stopping execution with a keyboard command\n');
   %fprintf('\nGive dbcont to continue the run, or dbquit to quit\n');
   fprintf(' ');
   error('ALOCDIRECT - Illegal input!!!')
end
useProb = 0;
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
% Set print levels
tomPrint        = PriLev - 5;
Prob.PriLevOpt  = 0;
ProbL.PriLevOpt = 0;
ProbL.WarmStart = 0;
if ~useProb
   N =ProbL.N;
   M = max(length(ProbL.c_L),length(ProbL.c_U))+ size(ProbL.A,1);
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
   probType = 3;
   N =ProbL.N;
   M = max(length(Prob.c_L),length(Prob.c_U))+ size(Prob.A,1);
   ProbL.optParam=optParamDef(localSolver,probType,N,N,M);
end
% optP=ProbL.optParam

if Prob.WarmStart == 1
   load(DIRECTFILE,'F');
   fPrev = length(F);
else
   fPrev = 0;
end

Prob.optParam.MaxFunc = maxFunc;
Prob.optParam.MaxIter = maxFunc;   % MaxIter is of less interest

% 1st step: Run glc with maxFunc evaluations
Prob.optParam.IterPrint = PriLev > 5;
r = tomRun(DIRECT,Prob,tomPrint);
Prob.optParam.IterPrint = 0;

load(DIRECTFILE,'feasible'); 

totfEval = r.FuncEv;
totcEval = r.ConstrEv;

% Avoid more testing of the structure
Prob.MENU = 1;

if ~feasible
   Prob.optParam.MaxFunc = maxFunc-1;
   ExitFlag = r.ExitFlag;
   while ~feasible & (totfEval < maxFunc3) & (ExitFlag == 0 | ExitFlag == 7)
      % Must try to get feasible
      Prob.WarmStart = 1;
      r = tomRun(DIRECT,Prob,tomPrint);
      totfEval = totfEval + r.FuncEv;
      totcEval = totcEval + r.ConstrEv;
      ExitFlag = r.ExitFlag;
      load(DIRECTFILE,'feasible'); 
   end
   Prob.WarmStart = 0;  % Reset WarmStart
end
if ~feasible & totfEval >= mFunc
   Result = r;
   Result.Solver          = 'ALOCDIRECT';
   Result.SolverAlgorithm = ...
    ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of major iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 7;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev) ...
                      ' f(x). No feasible point found'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('No feasible point found\n'); end
   return
end

Prob.WarmStart = 0;  % No WarmStart in the local solvers

% C and F could have been returned directly from glcDirect/glcSolve

if OLD
   load(DIRECTFILE,'C','glcfMin','F','ignoreIdx','fMinIdx','fMinEQ'); 
else
   load(DIRECTFILE,'C','glcfMin','F','ign','fMinIdx','fMinEQ'); 
   ignoreIdx = ign;
end
if PriLev > 0
   fprintf('ALOCDIRECT Phase 1: ')
   fprintf('Tried %4d. ',totfEval)
   if feasible
      fprintf('Found feasible best global point:',totfEval)
      fprintf('%20.10f\n',glcfMin)
      if fMinEQ > 0
         fprintf('                      ')
         fprintf('f(x) = %14.10f ',glcfMin-fMinEQ);
         fprintf('Deviation in equality constraints %12.10f\n',fMinEQ);
      end
   else
      fprintf('Found no feasible global point\n')
      if fMinEQ > 0
         fprintf('Deviation in equality constraints %12.10f\n',fMinEQ);
      end
   end
end
% Use fMinEQ nonzero to flag for local search to accept any feasible point
if ~feasible & fMinEQ == 0
   fMinEQ = 1;
end


% makeAllNorms returns the distance from any point
% to that points closest neighbour.
% minl1(i) = min (| C(:,i) - C(:,j)|, j=1:m, j~=i, m=size(C,2))

minl1=tomsol(29,C);
%minl1=makeAllNorms(C);

% 0.037 is an interesting number, often reappearing.

% meanDist=mean( minl1) is the average distance any point has to its closest neighbour
meanDist=mean(minl1);
% meanDist will be an upper bound on maxDist
% v is the set of points with that wont be clustered
vIdx=minl1>=meanDist;

v1=find(vIdx);
% v2 is the set of points being considered for clustering
v2=find(vIdx == 0);

[Fout Iout] = min(F(v1));

Cout = C(:,v1(Iout));

% Points that are to be clustered will no longer be considered for trisection
ignore            = zeros(length(F),1);
ignore(ignoreIdx) = 1;
ignore(v2)        = 1;
ignoreIdx         = find(ignore);

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

if length(v2) == 0
   Result = r;
   Result.Solver          = 'ALOCDIRECT';
   Result.SolverAlgorithm = ...
    ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of major iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 0;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev) ...
                      ' f(x). No points to cluster'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0 , fprintf('No points to cluster'); end
   return
end


v3  = cluster(C(:,v2),maxDist);
v4  = find(v3==0);
nPnt = 1+length(v4);

if PriLev > 2
   fprintf('                      ')
   fprintf('maxDist %8.5f ',maxDist)
   fprintf('nPnt %3d ',nPnt)
   fprintf('pnts Checked %d ',length(v2))
   fprintf('Total pnts %d ',size(C,2))
   fprintf('\n')
end

if nPnt > maxLocalTry
   while nPnt > maxLocalTry & maxDist < meanDist
      % Increase maxDist
      maxDist = incDist*maxDist;
      v3  = cluster(C(:,v2),maxDist);
      v4  = find(v3==0);
      nPnt = 1+length(v4);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(v2))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
elseif nPnt < minLocalTry
   maxDist = decDist*maxDist;
   while nPnt < minLocalTry & maxDist > maxDistMin
      v3S = v3;
      v4S = v4;
      v3  = cluster(C(:,v2),maxDist);
      v4  = find(v3==0);
      nPnt = 1+length(v4);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(v2))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
      maxDist = decDist*maxDist;
   end
   if nPnt > maxLocalTry
      % Too many clusters, use previous clustering
      v3      = v3S;
      v4      = v4S;
      maxDist = maxDist/decDist;
      nPnt    = 1+length(v4);
      if PriLev > 2
         fprintf('                      ')
         fprintf('Use previous clustering, nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(v2))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
   maxDist = maxDist/decDist;
elseif nPnt < maxLocalTry & incLocalTry
   maxDist = decDist*maxDist;
   while nPnt < maxLocalTry & maxDist > maxDistMin
      v3S = v3;
      v4S = v4;
      v3  = cluster(C(:,v2),maxDist);
      v4  = find(v3==0);
      nPnt = 1+length(v4);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(v2))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
      maxDist = decDist*maxDist;
   end
   if nPnt > maxLocalTry
      % Too many clusters, use previous clustering
      v3      = v3S;
      v4      = v4S;
      maxDist = maxDist/decDist;
      nPnt    = 1+length(v4);
      if PriLev > 2
         fprintf('                      ')
         fprintf('Use previous clustering, nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(v2))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
   maxDist = maxDist/decDist;
end

Fv  = zeros(nPnt,1);
c_k = zeros(n,nPnt);

if length(v4) >= 1
   ix = v2(v3(1:v4(1)-1));
   [Fv(1) idxMt] = min(F(ix));
   c_k(:,1)=C(:,ix(idxMt(1)));
   for i=2:length(v4)
       ix=v2(v3(v4(i-1)+1:v4(i)-1) );
       [Fv(i) idxMt] = min(F(ix));
       c_k(:,i)=C(:,ix(idxMt(1)));
   end
   ix=v2(v3(v4(length(v4))+1:length(v3)));
   [Fv(length(Fv)) idxMt] = min(F(ix));
   c_k(:,length(Fv))=C(:,ix(idxMt(1)));
else
   ix = v2(v3);
   [Fv(1) idxMt] = min(F(ix));
   c_k(:,1)=C(:,ix(idxMt(1)));
end
if Fout < max(Fv)
   % Add start from best non-clustered point
   Fv = [Fv;Fout*ones(size(Cout,2),1)];
   c_k = [c_k, Cout];
end

x_D = x_U - x_L;
x_k = tomsol(9,x_L,c_k,x_D);
%ixD = find(x_D > 0 & ~ixVars);
ixD = find(x_D > 0);

nPnt = size(c_k,2);
if PriLev == 2
   fprintf('                      ')
   fprintf('Number of clusters %d. maxDist %f\n',nPnt,maxDist)
end

if isempty(maxFunc2)
   left = mFunc - totfEval;
else
   left = mFunc - totfEval - maxFunc2;
end

if left < nPnt
   % Too few function values left, user must increase MaxFunc
   Result = r;
   Result.Solver          = 'ALOCDIRECT';
   Result.SolverAlgorithm = ...
    ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
   Result.Cluster.maxDist = maxDist;
   Result.Cluster.minDist = minl;
   Result.Cluster.minIdx  = v2;
   Result.Cluster.x_k     = x_k;
   Result.Cluster.f_k     = Fv;
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 1;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev)  ...
                      ' f(x). Max f(x) reached.'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('Max f(x) reached\n'); end
   return
end
if maxLocalTry <= 0
   % Only clustering should be made
   Result = r;
   Result.Solver          = 'ALOCDIRECT';
   Result.SolverAlgorithm = ['Constrained DIRECT with ' DIRECT '/ Cluster / NO Local Search'];
   Result.Cluster.maxDist = maxDist;
   Result.Cluster.minDist = minl;
   Result.Cluster.minIdx  = v2;
   Result.Cluster.x_k     = x_k;
   Result.Cluster.f_k     = Fv;
   Result.MinorIter= r.Iter;      % Number of iterations in glcDirect/glcSolve
   Result.Iter     = 1;           % Number of iterations
   Result.FuncEv   = totfEval;
   Result.ConstrEv = totcEval;
   Result.ExitFlag = 0;
   Result.ExitText = ['Tried ' num2str(totfEval+fPrev)  ...
                      ' f(x). Found ' num2str(nPnt) ' cluster points'];
   Result          = endSolve(Prob,Result);
   if PriLev > 0, fprintf('Clustering done, %d points found\n',nPnt); end
   return
end

%Prob.optParam.MaxIter=left/nPnt;
Prob.optParam.MaxIter=mIter;
Prob.optParam.MaxFunc=left/nPnt;

X = [];
Floc = [];
gradEv=0;
HessEv=0;
if useProb | 1
   % Take advantage of any derivatives in local search
   dc                     = Prob.FUNCS.dc;
   g                      = Prob.FUNCS.g;

   if isempty(g)
      if isempty(dc)
         DerLvl           = 0;
         cDiff            = 6;
         gDiff            = 6;
      else
         DerLvl           = 2;
         cDiff            = 0;
         gDiff            = 6;
      end
   elseif isempty(dc)
      DerLvl              = 1;
      cDiff               = 6;
      gDiff               = 0;
   else
      DerLvl              = 3;
      cDiff               = 0;
      gDiff               = 0;
   end
   
   switch lower(localSolver)
    case {'snopt','npsol','nlssol','minos'}
      % Use internal derivatives for the SOL solvers
      ProbL.NumDiff        = gDiff;
      ProbL.ConsDiff       = cDiff;
      ProbL.SOL.optPar(39) = DerLvl;
      %Prob.SOL.PrintFile=[localSolver '.txt'
      %Prob.SOL.SummFile= [localSolver 's.txt'
    otherwise
      ProbL.NumDiff        = gDiff > 0;
      ProbL.ConsDiff       = cDiff > 0;
   end
end
fMin = glcfMin;
xMin = r.x_k(:,1);
cDev = r.c_k;

if PriLev > 0
   fprintf('ALOCDIRECT Phase 2: ')
   fprintf('Number of initial points from clustering: %d\n',size(x_k,2))
end
   
ixBest = 0;
for i=1:size(x_k,2)
    ProbL.x_0 = x_k(:,i);
    if ~isempty(IntVars)
       x00 = max(x_L(IntVars),min(x_U(IntVars),ProbL.x_0(IntVars)));
       ProbL.x_L(IntVars)=x00;
       ProbL.x_U(IntVars)=x00;
    end
    if n == 1
       idx3 = find(C==c_k(:,i)*ones(1,size(C,2) )); 
    else
       idx3 = find(all( C==c_k(:,i)*ones(1,size(C,2)) )); 
    end
    % ProbL.SOL.PrintFile = ['snopt' num2str(i) '.txt'];
    rl   = tomRun(localSolver,ProbL,tomPrint);
    totfEval=totfEval + rl.FuncEv;
    totcEval=totcEval + rl.ConstrEv;
    gradEv=gradEv + rl.GradEv;
    HessEv=HessEv + rl.HessEv;
    if M > 0
       % Check if point feasible enough to be accepted
       h_L1 = consviolation(rl, PriLev-1);
       if h_L1 <= M*ProbL.optParam.cTol
          f_k  = rl.f_k;
          if PriLev > 1
             fprintf('                      ')
             fprintf('Local point accepted, f(x) = %20.10f ',f_k);
             fprintf('Local search %d\n',i);
          end
       else
          f_k = Inf;
       end
    else
       f_k = rl.f_k;
    end
    if (fMinEQ > 0 & ~isinf(f_k)) | f_k < fMin
       fMin = f_k;
       xMin = rl.x_k;
       cDev = rl.c_k;
       fMinEQ     = 0; % set better fMinEQ later on
       CMin = ProbL.x_L;
       CMin(ixD) = (rl.x_k(ixD) - x_L(ixD))./x_D(ixD);
       if ~isempty(IntVars)
          IMin = ProbL.x_L(IntVars);
       end
       ixMin = idx3(1);
       ixBest = i;
       if PriLev > 1 
          fprintf('                      ')
          fprintf('Step %d. ',i)
          fprintf('Improved best point with %s\n',localSolver)
          fprintf('                      ')
          fprintf('Old Best %18.14f. New Best %18.14f',glcfMin, f_k)
          fprintf('\n')
       end
    end
    if ~isempty(IntVars)
       ProbL.x_L(IntVars)=x_L(IntVars);
       ProbL.x_U(IntVars)=x_U(IntVars);
    end
    if PriLev > 5 
          fprintf('                      ')
          fprintf('End of local search # %d out of %d\n',i,size(x_k,2)); 
          fprintf('                      ')
          fprintf('----------------------------\n\n'); 
    end
    if totfEval > mFunc, break; end
end

if ixBest > 0
   if PriLev > 0 
      fprintf('ALOCDIRECT Phase 3: ')
      fprintf('Improved best point with local solver %s ',localSolver)
      fprintf('in local search %d\n',ixBest)
      fprintf('                    ')
      fprintf('Local Best: %19.10f ',fMin)
      if glcfMin >= 1E300
         fprintf('Global Best: %20.10f\n',Inf)
      else
         fprintf('Global Best: %20.10f\n',glcfMin)
      end
   end
   glcfMin    = fMin;
   if 1 
      Prob.MIP.fIP = glcfMin;
      Prob.MIP.xIP = xMin;
   else
      % Avoid this old way of changing the mat file
      %fMinIdx    = ixMin; % glcDirect/glcSolve returns empty x_k then
      ixMin      = fMinIdx(1);
      %ixMin      = 1;
      F(ixMin)   = fMin;
      C(:,ixMin) = CMin;
      if ~isempty(IntVars)
         load(DIRECTFILE,'iL','iU'); 
         iL(:,ixMin) = IMin;
         iU(:,ixMin) = IMin;
         save(DIRECTFILE,'iL','iU','-APPEND');
      end
      fMinIdx    = ixMin;
      fMinEQ     = 0; % set better fMinEQ later on
      Result.c_k = c_k;
      if OLD
         save(DIRECTFILE,'C','glcfMin','F','ignoreIdx','fMinIdx','fMinEQ','-APPEND')
      else
         ign = ignoreIdx;
         save(DIRECTFILE,'C','glcfMin','F','ign','fMinIdx','fMinEQ','-APPEND')
      end
      if feasible == 0
         if PriLev > 0 
            fprintf('Change DIRECT sample point %d to feasible ',ixMin);
            fprintf('using local minimum %d\n',ixBest);
         end
         feasible = 1;
         load(DIRECTFILE,'G')
         G(:,fMinIdx) = 0;
         save(DIRECTFILE,'feasible','G','-APPEND')
      end
   end
else
   fMin       = glcfMin;
end

Cd=C;
if isempty(maxFunc2)
   maxFunc = min(maxFunc1,max(n, mFunc - totfEval - 5*n));
else
   maxFunc = maxFunc2;
end
Prob.optParam.MaxFunc = maxFunc;
Prob.optParam.MaxIter = maxFunc;
%Prob.optParam.MaxFunc = 2000;
%Prob.optParam.MaxIter = 2000;

if length(ignoreIdx) == length(F)
   Phase4 = 0;
else
   Phase4 = 1;
   Prob.WarmStart = 1;
   % tomPrint=2
   % Prob.optParam.IterPrint = 1;
   % Prob.PriLevOpt = 2;
   Prob.optParam.IterPrint = PriLev > 5;
   Result = tomRun(DIRECT,Prob,tomPrint);

   %load(DIRECTFILE,'C');
   %C=setdiff(C',Cd','rows')';
   %plot(C(1,:),C(2,:),'m+')

   totfEval=totfEval + Result.FuncEv;
   totcEval=totcEval + Result.ConstrEv;
end

% *******************************************************************
if ~Phase4
   Phase5 = 0;
elseif isempty(maxFunc2)
   Phase5 = (mFunc - totfEval) > 0;
else
   Phase5 = maxFunc2 > 0;
end

load(DIRECTFILE,'feasible'); 
REJECT = 0;
if M > 0 & feasible & Phase5
   if isempty(Result.x_k)
      % No solution found
      Phase5 = 0;
   else
      Result.c_k = nlp_c(Result.x_k(:,1),Prob);
      % Check if point feasible enough to be accepted
      h_L1 = consviolation(Result, PriLev-1);
      if h_L1 > 10*M*ProbL.optParam.cTol
         % Infeasible second glcDirect/glcSolve solution
         if PriLev > 0 
            fprintf('                    ');
            fprintf('New point rejected due to infeasibility %f\n',h_L1);
         end
         Phase5 = 0;
         REJECT = 1;
      end
   end
end
if Phase5
   if glcfMin <= Result.f_k + fMinEQ
      % No improvement by glcDirect/glcSolve
      if PriLev > 0 
         fprintf('ALOCDIRECT-Phase 4: ')
         fprintf('Tried %4d. ',Result.FuncEv)
         fprintf('No improvement - Stop\n')
      end
      Phase5 = 0;
   else
      cDev   = Result.c_k;
   end
elseif Phase4 & feasible
   if REJECT
      if PriLev > 0 
         fprintf('ALOCDIRECT-Phase 4: ')
         fprintf('Tried %4d. ',Result.FuncEv)
         fprintf('Best point improved by global search ')
         fprintf('%16.10f (infeasible)\n',Result.f_k)
      end
   else
      if Result.f_k + fMinEQ < glcfMin
         glcfMin = Result.f_k + fMinEQ;
         fMin    = glcfMin;
         x_k     = Result.x_k;
         xMin    = x_k(:,1);
         cDev    = Result.c_k;
         if PriLev > 0 
            fprintf('ALOCDIRECT-Phase 4: ')
            fprintf('Tried %4d. ',Result.FuncEv)
            fprintf('Best point improved by global search ')
            fprintf('%20.10f\n',glcfMin)
         end
      end
   end
end
     
if Phase5
   % glcDirect/glcSolve improved the best point
   glcfMin = Result.f_k + fMinEQ;
   x_k     = Result.x_k;
   fMin    = glcfMin;
   xMin    = x_k(:,1);
   if PriLev > 0 
      fprintf('ALOCDIRECT Phase 4: ')
      fprintf('Tried %4d. ',Result.FuncEv)
      fprintf('Best point improved by global search ')
      fprintf('%20.10f\n',glcfMin)
      fprintf('Try local search from %d point(s)\n',size(x_k,2))
   end
   ProbL.optParam.MaxFunc = max(5*ProbL.N,mFunc-totfEval);
   %ProbL.optParam.MaxIter = ProbL.optParam.maxFunc;
   ProbL.optParam.MaxIter = mIter;
   ProbL.WarmStart = 0;
   ixBest = [];
   for i=1:min(100,size(x_k,2))
       ProbL.x_0 = x_k(:,i);
       if ~isempty(IntVars)
          x00 = max(x_L(IntVars),min(x_U(IntVars),ProbL.x_0(IntVars)));
          ProbL.x_L(IntVars)=x00;
          ProbL.x_U(IntVars)=x00;
       end
       rl   = tomRun(localSolver,ProbL,tomPrint+3);
       totfEval=totfEval + rl.FuncEv;
       totcEval=totcEval + rl.ConstrEv;
       gradEv=gradEv + rl.GradEv;
       HessEv=HessEv + rl.HessEv;
       if M > 0
          % Check if point feasible enough to be accepted
          h_L1 = consviolation(rl, PriLev-1);
          if h_L1 <= M*ProbL.optParam.cTol
             f_k = rl.f_k;
             c_k = rl.c_k;
             if PriLev > 1
                fprintf('                      ')
                fprintf('Local point accepted, f(x) = %20.10f\n',f_k);
             end
          else
             f_k = Inf;
          end
       else
          f_k = rl.f_k;
       end
       if f_k < fMin
          % Local point does improve
          fMin = f_k;
          CMin = ProbL.x_L;
          xMin = rl.x_k;
          cDev = rl.c_k;
          CMin(ixD) = (xMin(ixD) - x_L(ixD))./x_D(ixD);
          if ~isempty(IntVars)
             IMin = ProbL.x_L(IntVars);
          end
          ixMin = idx3(1);
          ixBest = i;
          if PriLev > 1 
             fprintf('                      ')
             fprintf('Step %d. ',i)
             fprintf('Improved best point with %s\n',localSolver)
             fprintf('                      ')
             fprintf('Old Best %18.14f. New Best %18.14f',glcfMin, f_k)
             fprintf('\n')
          end
       elseif f_k == fMin & ~isempty(ixBest)
          xMin = [xMin,rl.x_k];
          ixBest = [ixBest,i];
       end
       if ~isempty(IntVars)
          ProbL.x_L(IntVars)=x_L(IntVars);
          ProbL.x_U(IntVars)=x_U(IntVars);
       end
       if totfEval > mFunc, break; end
   end
else
   % Use Phase3 or Phase4 solution as best found
   Result.x_k = xMin;
   Result.f_k = fMin;
end
if Phase5 & fMin < glcfMin
   if PriLev > 0 
      fprintf('ALOCDIRECT Phase 5: ')
      fprintf('Improved best point with local solver %s ',localSolver)
      fprintf('in local search: ')
      fprintf(' %d',ixBest)
      fprintf('\n')
      fprintf('                      ')
      fprintf('Local Best: %20.10f ',fMin)
      if glcfMin >= 1E300
         fprintf('Global Best: %20.10f\n',Inf)
      else
         fprintf('Global Best: %20.10f\n',glcfMin)
      end
   end
   feasible   = 1;
   glcfMin    = fMin;
   if 1 
      Prob.MIP.fIP = glcfMin;
      Prob.MIP.xIP = xMin;
   else
      % Avoid this old way of changing the mat file
      % fMinIdx    = ixMin; 
      ixMin      = fMinIdx;
      F(ixMin)   = fMin;
      C(:,ixMin) = CMin;
      if OLD
         save(DIRECTFILE,'C','glcfMin','F','ignoreIdx','fMinIdx','-APPEND')
      else
         ign = ignoreIdx;
         save(DIRECTFILE,'C','glcfMin','F','ign','fMinIdx','-APPEND')
      end
      if ~isempty(IntVars)
         load(DIRECTFILE,'iL','iU'); 
         iL(:,ixMin) = IMin;
         iU(:,ixMin) = IMin;
         save(DIRECTFILE,'iL','iU','-APPEND');
      end
   end
   Result.x_k = xMin;
   Result.c_k = c_k;
elseif Phase5
   Result.x_k = x_k;
end

if glcfMin >= 1E300 & fMinIdx > 0
   % Safe guard if no feasible point found
   load(DIRECTFILE,'F');
   glcfMin = F(fMinIdx);
end

% ********************************************************************

Result.GradEv   = gradEv;
Result.HessEv   = HessEv;
Result.f_k      = glcfMin;     % Best function value
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
Result.Solver          = 'ALOCDIRECT';
Result.SolverAlgorithm = ...
 ['Constrained DIRECT with ' DIRECT '-Cluster-Local Search with ' localSolver];
Result.Cluster.maxDist = maxDist;
Result.Cluster.minDist = minl;
Result.Cluster.minIdx  = v2;
Result.Cluster.x_k     = x_k;
Result.Cluster.f_k     = Fv;
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
   
% MODIFICATION LOG:
%
% 060328  hkh  Written from glcCluster
% 060814  med  FUNCS used for callbacks instead
% 070222  hkh  Revise IntVars handling, use new format
