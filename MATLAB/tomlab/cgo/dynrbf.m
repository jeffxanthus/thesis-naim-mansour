% DYNRBF implements Dynamic RBF algorithms developed by Prof.Dr. K. Holmstrom
%
% DYNRBF is based on RBF algorithms presented in 
% 1. Bjorkman, Holmstrom: Global Optimization of Costly Nonconvext Functions 
%    Using Radial Basis Functions, Optimization and Engineering 1,373-397, 2000
% 2. Hans-Martin Gutmann: A radial basis function method for global 
%    optimization, Journal of Global Optimization, 19,201:207, 2001.
% 3. Hans-Martin Gutmann: Radial Basis Function Methods for Global Optimization,
%    Ph.D. Thesis, Cambridge University, Cambridge, UK, September 2001.
% and ...
%
% DYNRBF solves problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Any of the x could be set as integer valued
%
% f(x) are assumed to be a costly function
% c(x) are assumed to be cheaply computed
% If some subset c_I(x) are very costly, create f(x) as a penalty function as
%
%     NEW f(x) = f(x) + beta' * c_I(x), beta positive penalties
%
% Calling syntax:
%
% function Result = DYNRBF(Prob, varargin)
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say RBFF
%   FUNCS.c   The routine to compute the nonlinear constraint, say RBFC
%             A call to tomFiles.m or glcAssign.m sets these fields. 
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   WarmStart If true, >0, DYNRBF reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   PriLevOpt Print Level 
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%   PriLevSub Print Level in subproblem solvers, see help in snSolve, gnSolve 
%   user      User field used to send information to low-level functions
%   f_Low     Lower bound on the optimal function value. If defined, used to
%             restrict the target values into interval [f_Low,min(surface)]
% --------------------------------------------
% optParam    Structure in Prob, Prob.optParam 
% ---------------------------------------------
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1. See OUTPUT PRINTING below
%  MaxFunc    Maximal number of costly function evaluations. Default 300
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then MaxFunc = MaxFunc + nFunc
%  MaxIter    Maximal number of iterations used in the local optimization on
%             the response surface in each step. Default 1000, except for
%             pure IP problems, then max(GO.MaxFunc, MaxIter);
%  fGoal      Goal for function value, if empty or Inf not used
%  eps_f      Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%             See the output field maxTri.
%  bTol       Linear constraint tolerance
%  cTol       Nonlinear constraint tolerance
%
% ------------------
% Fields in Prob.CGO
% ------------------
%
% nSample      Number of sample points to be used in initial experimental 
%              design. The method is determined by the value of Percent.
%              nSample is used in different way dependent on Percent:
%              Percent > 100. 
%               nSample =[] gives nSample = (d+1)*(d+2)/2
%               nSample  >0 gives nSample = max(d+1,nSample)
%               nSample =0  gives nSample = d+1 (minimal number)
%               nSample  <0 gives the number of function evaluations
%                 as the result of doing abs(nSample) DIRECT iterations
%              Percent == 100. Field not used, nSample = size(Prob.CGO.X,2);
%              Percent > 0 (less than 100):
%               nSample =[] gives nSample = d+1 (minimal number)
%               nSample <=0 gives nSample = (d+1)*(d+2)/2
%               nSample  >0 gives nSample = max(d+1,nSample)
%              Percent == 0: 
%               nSample =[] gives all 2^d corner points 
%               nSample <=0 gives all 2^d corner points 
%               nSample  >0 nSample = max(d+1,nSample) of the corner points
%                           are randomly selected
%              Percent < 0: 
%               nSample =[] k = abs(percent) is used to determine # of points
%                           k      :  2  3  4  5  6  >6
%                           Points :  21 33 41 51 65 65
%               nSample <0  gives nSample = (d+1)*(d+2)/2
%               nSample ==0 gives nSample = d+1
%               nSample  >0 gives nSample = max(d+1,nSample)
%         
%              Otherwise nSample is predefined as:
%               Percent = -997: 2*d+2 + 1 (if Prob.AddMP)
%               Percent = -998:   d+1 + 1 (if Prob.AddMP)
%               Percent = -999:   d+1 + 1 (if Prob.AddMP)
%
% Percent      Type of strategy to get the initial sampled values:
%              Default Percent=-997 (if box-bounded problem)
%              Default Percent=-5000 (if constrained problem)
%
%              Percent > 100. Use determinstic global optimization methods for
%                the initial design. Current methods (all DIRECT methods)
%                101 = glcDirect
%                102 = glbDirect
%                103 = glcSolve
%                104 = glbSolve
%                105 = glcFast
%                106 = glbFast
%               nSample  >0 gives nSample = max(d+1,nSample)
%               nSample  <0 gives the number of function evaluations
%                 as the result of doing abs(nSample) DIRECT iterations
%              Percent == 100: User given initial points x, as a matrix of
%              points in CGO.X. Each column is one sampled point. 
%              If d = length(Prob.x_L), then size(X,1) = d, size(X,2) >= d+1
%              CGO.F should be defined as empty, or contain a vector of 
%              corresponding f(x) values. Any CGO.F value set as NaN will be 
%              computed by arbfmip.
%
%              Percent > 0 (less than 100): Random strategy, the Percent value 
%              gives the percentage size of an ellipsoid around the so far 
%              sampled points that the new points are not allowed in.
%              Range 1%-50%. Recommended values 10% - 20%.
%
%              Percent == 0: Initial points is the corner points of the box
%              Generates too many points if the dimension is high.
%
%              Percent < 0: Latin hypercube space-filling design. The value
%              of abs(Percent) should in principle be the dimension. The call
%              made is X = daceInit(round(abs(Percent)),Prob.x_L,Prob.x_U);
%              Let k = abs(Percent), then the number of points are:
%              k      :  2  3  4  5  6  >6
%              Points :  21 33 41 51 65 65
%              A more efficient number of points could be to instead use:
%              k = 1  :  Use (d+1)*(d+2)/2 points (enough to fit a quadratic)
%              d      : 1  2  3  4  5  6  7  8  9 10
%              Points : 3  6 10 15 21 28 36 45 55 66
%
%              Percent < -1000: Latin hypercube space-filling design, only keep
%              the points that fulfill the linear and nonlinear constraints. 
%              The algorithm will try up to M = abs(Percent)-1000 points,
%              stopping when it has got length(x_L)+1 feasible points
%
%              Percent == -997: Initial points are 2*d + 2 corners: 
%              the lower left corner x_L and
%              its d adjacent corners x_L + (x_U(i)-x_L(i))*e_i, i=1,...,d
%              and the upper right corner x_U and 
%              its d adjacent corners x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%              Percent == -998: Initial points are the upper right corner x_U
%              and its d adjacent corners x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%              Percent == -999: Initial points are the lower left corner x_L
%              and its d adjacent corners x_L + (x_U(i)-x_L(i))*e_i, i=1,...,d
%
% X            If Percent >= 100, a matrix of initial x values
%              One column for every x value. size(X,2) >= dim(x)+1 needed
% F            If Percent >= 100, a vector of initial f(x) values.
%              If any element is set to NaN, DYNRBF will compute f(x)
% CX           If Percent >= 100, optionally a matrix of nonlinear 
%              constraint c(x) values. If nonempty, then 
%              size(CX,2) == size(X,2). If any element in a column i is set as
%              NaN, the vector c(x) = CX(:,i) will be recomputed
%
% RandState    If >=0, rand('state',RandState) is set to initialize the
%              pseudo-random generator
%              if < 0, rand('state',sum(100*clock)) is set to give a new set
%              of random values each run
%              Default RandState = 0
%              If isnan(RandState), the random state is not initialized 
%              RandState will influence if a stochastic initial experimental
%              design is applied, see input Percent and nSample.
%              RandState will also influence if using the multiMin solver, but
%              the random state seed is not reset in multiMin.
%              The state of the random generator is saved in rngState, and the
%              random generator is reinitialized with this state if WarmStart
%              
%
% idea         Type of search strategy on the surface, idea 1 (cycle of N+1=5
%              points as default, with different target values fnStar) or 
%              idea 2 (cycle of 4 points in alpha, which implicitly gives 
%              the target value fnStar).
%              Default is idea 1, cycle length N = Prob.CGO.N.
% rbfType      Type of radial basis function.
%              1-Thin Plate Spline, 2-Cubic (Default).
% SCALE        0-Original search space (Default if any integer values)
%              1-Transform search space to unit cube (Default if no integers).
% PLOT         0-No plotting (Default), 1-Plotting sampled points.
% REPLACE      0-No replacement (Default for constrained problems)
%              1-Large function values are replaced by the median
%              >1 - Large values Z are replaced by new values
%              Replacement: Z:= FMAX + log10(Z-FMAX), where
%              FMAX = 10^REPLACE, if min(F) < 0 
%              FMAX = 10^(ceil(log10(min(F)))+REPLACE), if min(F) >= 0
%              (Default REPLACE = 5 if no linear or nonlinear constraints)
%              New replacement in every iteration, because min(F) may change
% LOCAL        0-No local searches after global search
%              If RBF surface is inaccurate, might be an advantage
%              1-Local search from best points after global search. If best
%              function values are equal, up to 20 local searches are done.
%              LOCAL is always set to 0 for multiMin and OQNLP.
%
% N            Cycle length in idea 1 (Default N=4 for fStarRule 1 and 2,
%              Default N=1 for fStarRule 3) or in idea 2 (always N=3)
%
% infStep      If =1, add search step with target value -inf first in cycle
%              Default 0.
% AddMP        If = 1, add the midpoint as extra point in the corner strategy or
%              the adjacent corner strategies Percent=0,-997,-998,-999. 
%              Default 0, except if Percent is any of 0,-997,-998,-999.
% AddSurfMin   Add up to AddSurfMin interior local minima on RBF surface as 
%              search points, based on estimated Lipschitz constants 
%              AddSurfMin=0 implies no additional minimum added (Default).
%              Only possible if globalSolver = 'multiMin'
%              Test for additional minimum in local step (modN == N)
%              modN = -2,-3,-4,... are iteration steps with these search points
% TargetMin    Which minimum of several to pick in target value problem
%              =0 Use global minimum
%              =1 Use best interior local minima, if none use global minimum
%              =2 Use best interior local minima, if none use RBF interior min
%              =3 Use best minimum with lowest number of coefficients on bounds
%              Default TargetMin = 3
% fStarRule    Global-Local search strategy in idea 1. N = cycle length
%              Define min_sn as global minimum on surface. fStar Target value
%              1: fStar = min_sn - ((N-(n-nInit))/N)^2*Delta_n (Default)
%              2: fStar = min_sn -  (N-(n-nInit))/N   *Delta_n
%              Strategy 1 and 2 depends on Delta_n estimate (see DeltaRule)
%              If infStep true, addition of -inf-step first in cycle
%              3: fStar = -inf-step, min_sn-k*0.1*|min_sn| k=N,...,0
%
%              Strategy names in Gutmanns thesis: III, II, I
%
% eps_sn       Relative tolerance used to test if minimum of surface, min_sn, is 
%              sufficiently lower than the best point found (fMin).  Default = 1E-7
% DeltaRule    1 = Skip large f(x) when computing f(x) interval Delta
%              0 = Use all points. Default 1.
% globalSolver Solver used for global optimization on the RBF surface
%              If the globalSolver is glcCluster, the fields
%              Prob.GO.maxFunc1, Prob.GO.maxFunc2 and Prob.GO.maxFunc3 are used
%              See the help for maxFunc1, maxFunc2, maxFunc3 in glcCluster
% localSolver  Solver used for local optimization on the RBF surface
%
% MaxCycle     Max number of cycles without progress before stopping. Default 10
%
% varargin     Additional parameters to DYNRBF are sent to the costly f(x)
%
% ---------------------------------------------------------
% Fields in Prob.GO (Default values are set for all fields)
% ---------------------------------------------------------
%
% MaxFunc      Maximal number of function evaluations in each global search
% MaxIter      Maximal number of iterations in each global search
% maxFunc1     glcCluster parameter, max function evaluations 1st call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc2     glcCluster parameter, max function evaluations 2nd call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc3     glcCluster parameter, max sum of function evaluations in 
%              repeated 1st calls trying to get feasible
%              Only used if globalSolver is glcCluster, see help globalSolver
% localSolver  The local solver used by glcCluster
% DIRECT       DIRECT method used in glcCluster: glcSolve or glcDirect(Default)
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:  
%             If empty, all variables are assumed non-integer 
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% If MIP problems are solved then the only subsolvers working are glcCluster,
% and OQNLP (for both the local and global subproblem)
% e.g.
% Prob.CGO.globalSolver = 'oqnlp';
% Prob.CGO.localSolver  = 'oqnlp';
% will use the OQNLP solver, a license for Tomlab /OQNLP is needed
%
% For pure IP problems, only glcSolve and glcDirect (Default) may be used. Set
% Prob.CGO.globalSolver = 'glcSolve'; to use glcSolve, otherwise glcDirect is used
%
% ---------------------------------------------
% OUTPUT PARAMETERS
% ---------------------------------------------
%
% Result    Structure with results from optimization
%  x_k      Matrix with the best points as columns, f(x_k) == f_k.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag Always 0, except
%           1 = Initial interpolation failed, normally because too huge f(x)
%  Inform   0 = Normal termination
%           1 = Function value f(x) is less than fGoal
%           2 = Error in function value f(x), abs(f-fGoal) <= fTol, fGoal=0
%           3 = Relative Error in function value f(x) is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           4 = Failure in global sub solvers, same point sampled over and over
%               Try lower cTol and or bTol, or other subsolver
%           7 = All feasible integers tried
%           8 = No progress for MaxCycle*(N+1)+1 function evaluations 
%               (>MaxCycle cycles, input CGO.MaxCycle)
%           9 = Max CPU Time reached
%
% To make a warm start possible, DYNRBF saves the following information in
% the file cgoSave.mat:
% Name      Name of the problem
% O         Matrix with sampled points (in original space)
% X         Matrix with sampled points (in unit space if SCALE == 1)
% F         Vector with function values
% F_m       Vector with function values (replaced)
% nInit     Number of initial points >= d+1 (2^d if center points)
% Fpen      Vector with function values + additional penalty if infeasible
% fMinIdx   Index of the best point found
% rngState  Current state of the random number generator used
%
% ---------------------------------------------
% OUTPUT PRINTING (IterPrint == 1 or PriLev > 0)
% ---------------------------------------------
% --- Row 1 
% Iter      Number of iterations
% n         Number of trial x, n-Iter is number of points in initial design
% nFunc     Number of costly f(x) computed, nFunc <= n, n-nFunc = rejected pnts
% Cycle     Cycle steps global to local. infStep is marked -1, 
%           0 to N-1 are global steps. Last step in cycle, N, is local search
% fnStar    Target value fn_star (alpha for idea 2).
% fGoal     Goal value (if set)
% fMin      Best f(x) found so far. E.g. at 27/It 12 means n=27, Iter=12 
%           fMinI means the best f(x) is infeasible
%           fMinF means the best f(x) is feasible (also integer feasible) 
%           
% fNew      f(xNew), the costly function value for point tried in current Iter
% In iteration 0 (if global optimum known and given in Prob.x_opt):
% dXO       Minimal distance from global optimum to closest point of all 
%           sampled points X in experimental design
% SumXO     Sum of distances from global optimum to all 
%           sampled points X in experimental design
% doO       Distance from xBest with best f(x) in experimental design 
%           to global optimum
% --- Row 2 (All distances are in SCALED space [0,1]^d ) 
% min_sn    Minimum on RBF surface, obtained at point min_sn_y
% [.]       Number of variables x active on bound at min_sn_y
%           Distances from min_sn_y to:
% doX       Minimal distance to closest point of all previous sampled points X
% doM       Distance to current best point found xMin, f(xMin) = fMinF
% doO       Distance to global optimum (if Prob.x_opt specified)
% [.]       Number of variables x active on bound at new point xNew
%           Distances from xNew to:
% doX       Minimal distance to closest point of all previous sampled points X
% doM       Distance to current best point found xMin, f(xMin) = fMinF
% doO       Distance to global optimum (if Prob.x_opt specified)
%
% SoO       Distance from min_sn_y to new point xNew
% 
% --- Row 3 
% xNew      The point tried in the current iteration
%
% USAGE:
%
% Let the name of the problem be "RBFF Test"
% The function RBFF is best written as
%     function f = RBFF(x, Prob) 
% Then any information, say u and W is easily sent to RBFF (and the constraint
% function RBFC, if present) using the Prob structure. See the example below. 
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%    Prob   = glcAssign('RBFF',x_L,x_U,'RBFF Test',A,b_L,b_U,'RBFC',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown on the two next lines
%    Prob.optParam.MaxFunc = 300;  % Change max number of function evaluations 
%    Prob.optParam.MaxIter = 2000; % Change local iterations to 2000
%    Prob.GO.MaxFunc       = 1000; % 1000 global function values
%
% Direct solver call:
%    Result = DYNRBF(Prob);
%    PrintResult(Result);
%           
% Driver call, including printing with level 2:
%      Result = tomRun('DYNRBF',Prob,2);
%           
% The user function RBFF is written as
%
%    function f = RBFF(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    f = "some function of x, u and W"
%
% It is also possible to use the function format
%    function f = RBFF(x) 
% but then any additional parameters must be sent as global variables.
%
% The user function RBFC, computing the nonlinear constraints, is written as
%
%    function c = RBFC(x, Prob) 
%    u = Prob.user.u; W = Prob.user.W;
%    c = "some vector function of x, V and W"
%
% Note! If RBFF has the 2nd input argument Prob, also RBFC must have that.
%
% To make a restart, just set the restart flag, and call DYNRBF once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('DYNRBF',Prob,2);

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2007 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Oct 8, 2000. Last modified Oct 10, 2007.

function Result = dynrbf(Prob, varargin)

if nargin < 1
   error('DYNRBF needs one parameter, the structure Prob');
end
global NLP_x NLP_f NARG

time  = fix(clock);

DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt
SolveInf    = 0;  % If to solve the inf problem first in every step, time consuming

solvType=checkType('glc');

Prob=ProbCheck(Prob,'DYNRBF',solvType);

Prob = iniSolve(Prob,solvType,0,0);

Result                 = ResultDef(Prob);
Result.Solver          = 'DYNRBF';
Result.SolverAlgorithm = 'Radial Basis Function Interpolation';

% Pick up bounds from Prob structure and check if OK
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds

if isempty(x_L) | isempty(x_U)
   disp('DYNRBF requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'DYNRBF requires both lower and upper variable bounds';
   Result.DIGIT    = [];
   Result.CGO.Its  = [];
   Result=endSolve(Prob,Result);
   return;
end
if ~(all(isfinite(x_L)) & all(isfinite(Prob.x_U)))
   disp('DYNRBF only solves box bounded problems, where both');
   disp('lower bounds Prob.x_L and upper bounds Prob.x_U are finite');
   Result.ExitFlag = 1;
   Result.ExitText = 'Problem not box bounded, variable bounds on x not finite';
   Result.DIGIT    = [];
   Result.CGO.Its  = [];
   Result=endSolve(Prob,Result);
   return;
end

% Pick up input variables from Prob structure
MaxCPU    = Prob.MaxCPU;
PriLev    = Prob.PriLevOpt;          % Print level
f_Low     = Prob.f_Low;              % Lower bound on f
MaxIter   = Prob.optParam.MaxIter;   % Iterations in local suboptimization
MaxFunc   = Prob.optParam.MaxFunc;   % Number of costly function evaluations
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration

fTol      = Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
fGoal     = Prob.optParam.fGoal;     % Goal f(x) for the optimization
epsRank   = Prob.optParam.eps_Rank;  % Rank tolerance in qr-decomposition
bTol      = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance

epsX      = 1.7E-6;

if isfield(Prob,'PriLevSub') 
   PriSub = Prob.PriLevSub;
else
   PriSub = [];
end
if isempty(PriSub), PriSub = 0; end

%xTol      = Prob.optParam.eps_x;     % Tolerance for rectangle sizes
%                                     % (scaled to (0,1) )

%eps_sn    = 1E-6; % Tolerance for min_sn

x_D   = x_U - x_L;
d     = length(x_L);  % Problem dimension
nMax  = Inf;
if sum(x_D) == 0
   %error('All variables are fixed')
end

% Safeguard
if isempty(MaxFunc), MaxFunc = 300; end
if MaxFunc < 0
   MaxFunc = 300;
end
MaxFunc = min(MaxFunc,5000); % Safe guard to avoid crash
if isempty(MaxIter), MaxIter = 1000; end
if MaxIter < 0
   MaxIter = 1000;
end
if isempty(IterPrint), IterPrint = 1; end

if isempty(Prob.CGO)
   idea      = []; rbfType   = []; SCALE        = []; PLOT        = [];
   REPLACE   = []; Percent   = []; globalSolver = []; localSolver = [];
   LOCAL     = []; N         = []; infStep      = []; AddMP       = []; 
   fStarRule = []; DeltaRule = []; RandState    = []; MaxCycle    = [];
   nSample   = []; eps_sn    = []; AddSurfMin   = []; TargetMin   = []; 


else

   if isfield(Prob.CGO,'idea') 
      idea = Prob.CGO.idea;
   else
      idea = [];
   end
   if isfield(Prob.CGO,'rbfType')
      rbfType = Prob.CGO.rbfType;
   else
      rbfType = [];
   end
   if isfield(Prob.CGO,'SCALE')
      SCALE = Prob.CGO.SCALE;
   else
      SCALE = [];
   end
   if isfield(Prob.CGO,'PLOT')
      PLOT = Prob.CGO.PLOT;
   else
      PLOT = [];
   end
   if isfield(Prob.CGO,'REPLACE')
      REPLACE = Prob.CGO.REPLACE;
   else
      REPLACE = [];
   end
   if isfield(Prob.CGO,'LOCAL')
      LOCAL = Prob.CGO.LOCAL;
   else
      LOCAL = [];
   end
   if isfield(Prob.CGO,'N')
      N = Prob.CGO.N;
   else
      N = [];
   end
   if isfield(Prob.CGO,'infStep')
      infStep = Prob.CGO.infStep;
   else
      infStep = [];
   end
   if isfield(Prob.CGO,'AddMP')
      AddMP = Prob.CGO.AddMP;
   else
      AddMP = [];
   end
   if isfield(Prob.CGO,'AddSurfMin')
      AddSurfMin = Prob.CGO.AddSurfMin;
   else
      AddSurfMin = [];
   end
   if isfield(Prob.CGO,'TargetMin')
      TargetMin = Prob.CGO.TargetMin;
   else
      TargetMin = [];
   end
   if isfield(Prob.CGO,'fStarRule')
      fStarRule = Prob.CGO.fStarRule;
   else
      fStarRule = [];
   end
   if isfield(Prob.CGO,'eps_sn')
      eps_sn = Prob.CGO.eps_sn;
   else
      eps_sn = [];
   end
   if isfield(Prob.CGO,'DeltaRule')
      DeltaRule = Prob.CGO.DeltaRule;
   else
      DeltaRule = [];
   end
   if isfield(Prob.CGO,'Percent')
      Percent = Prob.CGO.Percent;
   else
      Percent = [];
   end
   if isfield(Prob.CGO,'nSample')
      nSample = Prob.CGO.nSample;
   else
      nSample = [];
   end
   if isfield(Prob.CGO,'RandState')
      RandState = Prob.CGO.RandState;
   else
      RandState = [];
   end
   if isfield(Prob.CGO,'globalSolver')
      globalSolver = deblank(Prob.CGO.globalSolver);
   else
      globalSolver = [];
   end
   if isfield(Prob.CGO,'localSolver')
      localSolver = deblank(Prob.CGO.localSolver);
   else
      localSolver = [];
   end
   if isfield(Prob.CGO,'MaxCycle') 
      MaxCycle = Prob.CGO.MaxCycle;
   else
      MaxCycle = [];
   end
end
if isempty(Prob.GO)
   GOMaxFunc = []; GOMaxIter     = []; maxFunc1 = [];  maxFunc2 = [];
   maxFunc3  = []; GOlocalSolver = []; DIRECT   = [];
else
   if isfield(Prob.GO,'MaxFunc')
      GOMaxFunc = Prob.GO.MaxFunc;
   else
      GOMaxFunc = [];
   end
   if isfield(Prob.GO,'MaxIter')
      GOMaxIter = Prob.GO.MaxIter;
   else
      GOMaxIter = [];
   end
   if isfield(Prob.GO,'maxFunc1')
      maxFunc1 = Prob.GO.maxFunc1;
   else
      maxFunc1 = [];
   end
   if isfield(Prob.GO,'maxFunc2')
      maxFunc2 = Prob.GO.maxFunc2;
   else
      maxFunc2 = [];
   end
   if isfield(Prob.GO,'maxFunc3')
      maxFunc3 = Prob.GO.maxFunc3;
   else
      maxFunc3 = [];
   end
   if isfield(Prob.GO,'localSolver')
      GOlocalSolver = Prob.GO.localSolver;
   else
      GOlocalSolver = [];
   end
   if isfield(Prob.GO,'DIRECT')
      DIRECT = Prob.GO.DIRECT;
   else
      DIRECT = [];
   end
end
% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(d,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > d)
      error('DYNRBF: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

if ~isempty(IntVars)
   if length(IntVars) == d
      % Pure IP problem
      nMax = prod(1+(x_U-x_L));
   else
      ix = ones(d,1);
      ix(IntVars) = 0;
      if all(x_L(find(ix))==x_U(find(ix)))
         % All continuous variables are fixed
         nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
      end
   end
   if length(IntVars) == d
      % Pure IP problem
      if strcmpi(globalSolver,'glcSolve')
         globalSolver = 'glcSolve';
         localSolver  = 'glcSolve';
      elseif strcmpi(globalSolver,'oqnlp')
         globalSolver = 'oqnlp';
         localSolver  = 'oqnlp';
      elseif strcmpi(globalSolver,'glcFast')
         globalSolver = 'glcFast';
         localSolver  = 'glcFast';
      else
         globalSolver = 'glcDirect';
         localSolver  = 'glcDirect';
      end
      LOCAL        = 0;
      if isempty(GOMaxFunc)
         GOMaxFunc    = 5000;
      end
      if isempty(GOMaxIter)
         GOMaxIter    = 5000;
      end
      if isempty(MaxIter)
         MaxIter      = max(GOMaxFunc,GOMaxIter);
      else
         MaxIter      = max(MaxIter,GOMaxFunc);
      end
   elseif strcmpi(globalSolver,'oqnlp')
      if ~strcmpi(localSolver,'oqnlp') & ~strcmpi(localSolver,'glcCluster')
         localSolver  = 'oqnlp';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcDirect')
      if isempty(localSolver)
         localSolver  = 'glcDirect';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcDirect')
      if isempty(localSolver)
         localSolver  = 'glcFast';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcSolve')
      if isempty(localSolver)
         localSolver  = 'glcSolve';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'minlpSolve')
      if isempty(localSolver)
         localSolver  = 'minlpSolve';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'minlpBB')
      if isempty(localSolver)
         localSolver  = 'minlpBB';
      end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcCluster')
      globalSolver = 'glcCluster';
      if isempty(GOlocalSolver)
         if isempty(localSolver)
            GOlocalSolver = GetSolver('con',1,0);
         else
            GOlocalSolver  = localSolver;
         end
      end
      localSolver  = 'glcCluster';
      LOCAL = 0;
   elseif strcmpi(globalSolver,'multiMin')
      globalSolver = 'multiMin';
      if isempty(GOlocalSolver)
         if isempty(localSolver)
            GOlocalSolver = GetSolver('con',1,0);
         else
            GOlocalSolver  = localSolver;
         end
      end
      localSolver  = 'multiMin';
      LOCAL = 0;
   else % Default MINLP subsolver glcCluster
      globalSolver = 'glcCluster';
      if isempty(GOlocalSolver)
         if isempty(localSolver)
            GOlocalSolver = GetSolver('con',1,0);
         else
            GOlocalSolver  = localSolver;
         end
      end
      localSolver  = 'glcCluster';
   end
   SCALE = 0;
   % Safe guard bounds onto integer values
   x_L(IntVars)   = ceil(x_L(IntVars)); 
   x_U(IntVars)   = floor(x_U(IntVars)); 
   x_D            = x_U - x_L;
end
Prob.MIP.IntVars = IntVars;

nFunc     = 0;
nCon      = 0;        % Count number of calls to constraint routine, nlp_c
convflag  = 0;
nmax      = [];       % Only used in idea 1
CX        = [];       % Initial set of nonlinear constraint values
DEBUG     = 0;


% Default strategies for RBF input
if isempty(RandState),    RandState = 0; end
if isempty(idea),         idea = 1; end
if isempty(rbfType),      rbfType = 2; end
if isempty(SCALE)
   if all(x_L==0) & all(x_U==1)
      SCALE = 0;
   else
      SCALE = 1;
   end
end
if isempty(PLOT),         PLOT = 0; end
if isempty(LOCAL),        LOCAL = 1; end
if isempty(infStep),      infStep = 0; end
if isempty(AddSurfMin),   AddSurfMin = 0; end
if isempty(TargetMin),    TargetMin = 3; end
if isempty(fStarRule),    fStarRule = 1; end
if isempty(eps_sn),       eps_sn    = 1E-7; end
if isempty(DeltaRule),    DeltaRule = 1; end
%if isempty(globalSolver) 
%   globalSolver = 'multiMin'; 
%   localSolver  = 'multiMin';
%end
% if isempty(globalSolver), globalSolver = 'glcDirect'; end
if isempty(globalSolver), globalSolver = 'multiMin'; end
if isempty(localSolver),  localSolver = GetSolver('con',1,0); end
if isempty(MaxCycle),     MaxCycle = 10; end

if idea == 1      % Cycle length for direct fnStar strategy
   if fStarRule == 3
      infStep = 1; % Must have infStep in strategy 3, (Gutmann I)
      if isempty(N), N = 1; end
   else
      if isempty(N), N = 4; end
   end
elseif idea == 2  % Cycle length for indirect alpha strategy
   %if isempty(N), N = 3; end
   N = 3; % Always N = 3
end


% Possibly amplify fDiff to generate wider range of target values
AmpFac = 1.0;


dLin = size(Prob.A,1);
if dLin > 0
   d1 = length(Prob.b_L);
   % Must adjust bounds to have equal length
   if d1 < dLin
      Prob.b_L = [Prob.b_L;-inf*ones(dLin-d1,1)];
   end
   d1 = length(Prob.b_U);
   if d1 < dLin
      Prob.b_U = [Prob.b_U;inf*ones(dLin-d1,1)];
   end
end
dCon = max(length(Prob.c_L),length(Prob.c_U));
if dCon > 0
   d1 = length(Prob.c_L);
   % Must adjust bounds to have equal length, to add extra EGO constraint
   if d1 < dCon
      Prob.c_L = [Prob.c_L;-inf*ones(dCon-d1,1)];
   end
   d1 = length(Prob.c_U);
   if d1 < dCon
      Prob.c_U = [Prob.c_U;inf*ones(dCon-d1,1)];
   end
end
% Default values for GO structure
if dLin + dCon > 0
   if isempty(REPLACE),      REPLACE = 0; end
   if isempty(Percent),      Percent = -5000; end
   if strcmpi(globalSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      if isempty(maxFunc1),     maxFunc1  = 1000; end
      if isempty(maxFunc2),     maxFunc2  = 2000; end
      if isempty(maxFunc3),     maxFunc3  = 3000; end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcDirect')
      if isempty(GOMaxFunc),    GOMaxFunc = max(3000,300*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,300*d); end
   elseif strcmpi(globalSolver,'glcFast')
      if isempty(GOMaxFunc),    GOMaxFunc = max(3000,300*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,300*d); end
   elseif strcmpi(globalSolver,'glcSolve')
      if isempty(GOMaxFunc),    GOMaxFunc = max(1500,150*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1500,150*d); end
   elseif strcmpi(globalSolver,'multiMin')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      LOCAL = 0;
   else
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1000,1000*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
   end
   if strcmpi(localSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      if isempty(maxFunc1),     maxFunc1  = 1000; end
      if isempty(maxFunc2),     maxFunc2  = 2000; end
      if isempty(maxFunc3),     maxFunc3  = 3000; end
      MaxIter = max(MaxIter,GOMaxIter);
   end
else
   if isempty(REPLACE),      REPLACE = 5; end
   if isempty(Percent),      Percent = -d; end
   if strcmpi(globalSolver,'glcCluster')
      if isempty(GOMaxFunc),    GOMaxFunc = max(5000,500*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,250*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
      % Settings in EGO. Better?
      %if isempty(GOMaxFunc),    GOMaxFunc = max(5000,1000*d); end
      %if isempty(GOMaxIter),    GOMaxIter = max(3000,500*d); end
      %if isempty(maxFunc1),     maxFunc1  = 200+d*200; end
      %if isempty(maxFunc2),     maxFunc2  = 0; end
      %if isempty(maxFunc3),     maxFunc3  = maxFunc1; end
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcDirect')
      if isempty(GOMaxFunc),    GOMaxFunc = max(2000,200*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(2000,200*d); end
   elseif strcmpi(globalSolver,'glcFast')
      if isempty(GOMaxFunc),    GOMaxFunc = max(2000,200*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(2000,200*d); end
   elseif strcmpi(globalSolver,'glcSolve')
      if isempty(GOMaxFunc),    GOMaxFunc = max(1000,100*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(1000,100*d); end
   elseif strcmpi(globalSolver,'multiMin')
      if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
      LOCAL = 0;
   else
      if isempty(GOMaxFunc),    GOMaxFunc = max(5000,500*d); end
      if isempty(GOMaxIter),    GOMaxIter = max(3000,250*d); end
      if isempty(maxFunc1),     maxFunc1  = 500; end
      if isempty(maxFunc2),     maxFunc2  = 500; end
      if isempty(maxFunc3),     maxFunc3  = 500; end
   end
end
if isempty(AddMP),        
   if any(Percent == [0,-997,-998,-999])
      % Add the midpoint for the corner and adjacent corner strategies
      AddMP = 1; 
   else
      AddMP = 0; 
   end
end
% NOTE! Hard coded constants determining current cycle strategy.
% Should be further tested and developed to see which one is best.
% The old DYNRBF
OldStrat = 0; % If false, use a new Idea 1 median strategy
if dLin + dCon > 0
   OldAlpha = 0; % If false, use a new Idea 2 alpha  strategy
else
   OldAlpha = 1; % If false, use a new Idea 2 alpha  strategy
end

% Scaling parameters
if SCALE
   LOWER = 0;
   UPPER = 1;
   x_LL = LOWER*ones(d,1);
   x_UU = UPPER*ones(d,1);
   x_DD = x_UU-x_LL;
else
   x_LL = x_L;
   x_UU = x_U;
   x_DD = x_D;
end

if Percent > 0 & Percent < 100
   if isempty(nSample)
      nSample = d+1;
   elseif nSample <= 0
      nSample = (d+1)*(d+2)/2;
   else
      nSample = max(d+1,nSample);
   end
elseif Percent == 0
   if isempty(nSample)
      %nSample = [];
   elseif nSample <= 0
      nSample = [];
   else
      nSample = max(d+1,nSample);
   end
elseif Percent < 0
   if isempty(nSample)
      %nSample = [];
   elseif nSample < 0
      nSample = (d+1)*(d+2)/2;
   elseif nSample == 0
      nSample = d+1;
   else
      nSample = max(d+1,nSample);
   end
end

% Set possibly redefined SCALE into Prob
Prob.CGO.SCALE     = SCALE;
Prob.CGO.RandState = RandState;

if checkLicense('oqnlp')
   if isempty(IntVars)
      backupSolver1      = 'multiMin';
      backupSolver2      = 'oqnlp';
   else
      backupSolver1      = 'oqnlp';
      backupSolver2      = 'multiMin';
   end
else
   backupSolver1      = 'multiMin';
   backupSolver2      = [];
end

% Save input parameters in output structure Result.CGO
Result.CGO = struct('RandState',RandState,...
      'Percent',Percent,'nSample',nSample,'fTol',fTol,...
      'infStep',infStep,'N',N,'AddMP',AddMP, 'idea',idea, 'rbfType',rbfType, ...
      'fStarRule',fStarRule, 'eps_sn', eps_sn,'DeltaRule',DeltaRule, ...
      'REPLACE',REPLACE,'SCALE',SCALE,'MaxCycle',MaxCycle,'LOCAL',LOCAL,...
      'AddSurfMin',AddSurfMin,'TargetMin',TargetMin,'f_Low',f_Low,...
      'localSolver',localSolver,'globalSolver',globalSolver,...
      'GOMaxFunc',GOMaxFunc,'GOMaxIter',GOMaxIter,...
      'maxFunc1',maxFunc1,'maxFunc2',maxFunc2,'maxFunc3',maxFunc3,...
      'backupSolver1',backupSolver1,'backupSolver2',backupSolver2);

% Use optimal value for printing, if given
if isempty(Prob.x_opt)
   x_opt = [];
else
   x_opt = Prob.x_opt(1,:);
end
if SCALE
   if ~isempty(x_opt)
      xOptS   = (x_opt(:)-x_L)./x_D; 
   else
      xOptS   = [];
   end
else
   xOptS   = x_opt(:);
end

%
%  STEP 1, Initialization
%

if Prob.WarmStart
   % Restart with values from previous run.

   WS = DefPar(Prob.CGO,'WarmStartInfo');
   if ~isempty(WS)
      Name     = WS.Name;
      O        = WS.O;
      F        = WS.F;
      X        = WS.X;
      F_m      = WS.F_m;
      nInit    = WS.nInit;
      Fpen     = WS.Fpen;
      fMinIdx  = WS.fMinIdx;
      rngState = WS.rngState;
   else
      load('cgoSave.mat','Name','O','F','X','F_m','nInit','Fpen','fMinIdx',...
           'rngState');
   end
   rand('state',rngState);

   Name1 = deblank(Prob.Name);   % Name for the problem

   if strcmp(Name1,Name)
      [d n] = size(X);
      nFunc = n; % Total count of function evaluations
   
      if PriLev > 0
         fprintf('\n Restarting with %d sampled points from previous run\n',n);
      end
   else
      nFunc = 0; % Total count of function evaluations
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files cgoSave.mat?\n');
      end
   end

   % TOMSOL INIT, send F_m to Fortran, not F
   %clear tomsol;
   clear('tomsol');
   control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
   if control < 0
       fprintf('Initial interpolation failed');
       tomsol(25); % Deallocates memory
       Result.ExitFlag = 1;
       Result.ExitText = 'Initial interpolation failed';
       Result.DIGIT    = [];
       Result.CGO.Its  = [];
       Result          = endSolve(Prob,Result);
       return      %Something is really wrong
   end
   
   % Add extra function evaluations 
   if MaxFunc <= nFunc
      MaxFunc = MaxFunc + nFunc;
   end
   
end

if ~Prob.WarmStart
   Name = deblank(Prob.Name);  % Problem name
      
   % Set pseudo random generator
   if isnan (RandState)
      % Do no initialization
   elseif length(RandState) >1 | RandState >= 0
      rand('state',RandState);
   else
      rand('state',sum(100*clock));
   end
   % Three DIFFERENT INITIALIZATION STRATEGIES, plus user defined
   [X,O,F,Fpen,C,nCon,nSample] = expDesign(nSample, Percent, AddMP, ...
          SCALE,NaN,[PriLev >0 | IterPrint], Prob, varargin{:});

   Result.CGO.nSample = nSample;

   % Set dimension d and number of points n
   [d n] = size(X);
       
   % Set initial number of points nInit to n
   nInit = n;

   if REPLACE > 1
      % Replace very large function values by log
      FMIN = min(F);
      if FMIN <= 0
         FMAX = 10^REPLACE;
      else
         FMAX = 10^(ceil(log10(FMIN))+REPLACE);
      end
      ix = find(F>FMAX);
      F_m = F;
      if ~isempty(ix)
         fprintf('REPLACE %d large values with ',length(ix));
         fprintf('FMAX+log10(F-FMAX+1). FMAX=%e\n',FMAX);
         F_m(ix) = FMAX + log10(F(ix)-FMAX+1);
      end
   elseif REPLACE == 1
      % Replace large function values by median(F)
      F_m = min(median(F),F);
   else
      F_m = F;
   end
 
   % TOMSOL INIT, send F_m to Fortran, not F
   %clear tomsol;
   clear('tomsol');
   control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
   if control < 0
       fprintf('Initial interpolation failed\n');
       tomsol(25) % Deallocates memory
       if REPLACE == 0
          REPLACE = 1;
          fprintf('Retry with REPLACE=1, i.e. using median(F) avoiding');
          fprintf(' very large f(x) values\n');
          F_m = min(median(F),F);
          control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
       end
       if control < 0
          fprintf('Initial interpolation failed with REPLACE=1\n');
          fprintf('Lowest function value %30.20f\n',min(F));
          fprintf('Median function value %30.20f\n',median(F));
          fprintf('When quote median(F) / min(F) too high, try scaling ')
          fprintf('or transforming F, e.g. log(F) if f(x) > 0\n');
          for i =1:length(F)
              fprintf('Function value %d is %30.20f\n',i,F(i));
          end
          Result.ExitFlag = 1;
          Result.ExitText = 'Initial interpolation failed with REPLACE=1';
          Result.CGO.Its  = [];
          Result.DIGIT    = [];
          Result          = endSolve(Prob,Result);
          return     %Something is really wrong
       end
   end
   nFunc = n; % Number of function evaluations.
end

dc = Prob.FUNCS.dc;

if SCALE > 0 & (dCon > 0 | dLin > 0)
   if dLin > 0
      A   = Prob.A*diag(x_D);
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
   else
      A   = [];
      b_L = [];
      b_U = [];
   end
   if dCon > 0
      snProb             = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [],[],[], A,b_L,b_U, 'rbf_c','rbf_dc','rbf_d2c',...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
      snProb.xD          = x_D;
      snProb.xL          = x_L;
      snProb.cNargin     = xnargin(Prob.FUNCS.c);
      if isempty(dc)
         snProb.dcNargin = 0;
      else
         snProb.dcNargin = xnargin(dc);
      end
      snProb.c           = Prob.FUNCS.c;
      snProb.dc          = Prob.FUNCS.dc;
      snProb.d2c         = Prob.FUNCS.d2c;
   else
      snProb             = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [], [], [], A, b_L, b_U); 
   end
   snProb.SCALE          = SCALE;

elseif dCon == 0 & dLin == 0
   snProb                = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn'); 
else
   snProb                = conAssign('sn_f','sn_g',[],[],x_LL,x_UU,'RBFsn', ...
                           [], [], [], Prob.A, Prob.b_L, Prob.b_U, ...
                           Prob.FUNCS.c, Prob.FUNCS.dc, Prob.FUNCS.d2c, ...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
end

snProb.dDim               = d;
optParam                  = optParamDef(localSolver,solvType,d,dCon,dCon+dLin);
snProb.optParam           = optParam;
snProb.optParam.MaxIter   = MaxIter;
if strcmpi(localSolver,'glcCluster')
   snProb.optParam.MaxFunc   = GOMaxFunc;
else
   snProb.optParam.MaxFunc   = MaxIter*max(d,10)/10;  
end
snProb.optParam.IterPrint = DebugPriLev > 0;
snProb.optParam.cTol      = cTol;
snProb.optParam.bTol      = bTol;
snProb.GradTolg           = Prob.GradTolg;
snProb.GradTolH           = Prob.GradTolH;
snProb.GradTolJ           = Prob.GradTolJ;
snProb.SOL                = Prob.SOL;
snProb.GO                 = Prob.GO;
if isfield(Prob,'OQNLP')
   snProb.OQNLP           = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   snProb.KNITRO          = Prob.KNITRO;
end

snProb.PriLevOpt          = DebugPriLev;

% Send user info
if isfield(Prob,'user')
   snProb.user            = Prob.user;
end
snProb.MIP                = Prob.MIP;
snProb.uP                 = Prob.uP;
snProb.P                  = Prob.P;
snProb.mLin               = dLin;
snProb.mNonLin            = dCon;


if SCALE > 0 & (dCon > 0 | dLin > 0)
   if dLin > 0
      A   = Prob.A*diag(x_D);
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
   else
      A   = [];
      b_L = [];
      b_U = [];
   end
   if dCon > 0
      gnProb             = glcAssign('gn_f',x_LL,x_UU,'RBFgn',A,...
                           b_L,b_U,'rbf_c',Prob.c_L,Prob.c_U);
      gnProb.xD          = x_D;
      gnProb.xL          = x_L;
      gnProb.FUNCS.dc     = 'rbf_dc';
      gnProb.FUNCS.d2c    = 'rbf_d2c';
      gnProb.c           = Prob.FUNCS.c;
      gnProb.dc          = Prob.FUNCS.dc;
      gnProb.d2c         = Prob.FUNCS.d2c;
      gnProb.cNargin     = xnargin(Prob.FUNCS.c);
      if isempty(dc)
         gnProb.dcNargin = 0;
      else
         gnProb.dcNargin = xnargin(dc);
      end
      gnProb.SCALE       = SCALE;
   else
      gnProb             = glcAssign('gn_f',x_LL,x_UU,'RBFgn',A, b_L,b_U);
   end

elseif dCon == 0 & dLin == 0
   gnProb                = glcAssign('gn_f',x_LL,x_UU,'RBFgn');
else
   gnProb                = glcAssign('gn_f',x_LL,x_UU,'RBFgn',Prob.A,...
                           Prob.b_L,Prob.b_U,Prob.FUNCS.c,Prob.c_L,Prob.c_U);
   gnProb.FUNCS.dc        = Prob.FUNCS.dc;
end

% Send user info
if isfield(Prob,'user')
   gnProb.user           = Prob.user;
end
gnProb.MIP               = Prob.MIP;
gnProb.uP                = Prob.uP;
gnProb.P                 = Prob.P;
gnProb.dDim              = d;
gnProb.mLin              = dLin;
gnProb.mNonLin           = dCon;

optParam                 = optParamDef(globalSolver,solvType,d,dCon,dCon+dLin);
gnProb.optParam          = optParam;
gnProb.optParam.MaxIter  = GOMaxIter; 
gnProb.optParam.MaxFunc  = GOMaxFunc;  
gnProb.optParam.eps_Rank = epsRank;
gnProb.optParam.IterPrint = DebugPriLev > 0;
gnProb.optParam.cTol     = cTol;
gnProb.optParam.bTol     = bTol;
gnProb.GradTolg          = Prob.GradTolg;
gnProb.GradTolH          = Prob.GradTolH;
gnProb.GradTolJ          = Prob.GradTolJ;
gnProb.GO                = Prob.GO;
gnProb.SOL               = Prob.SOL;
if isfield(Prob,'OQNLP')
   gnProb.OQNLP          = Prob.OQNLP;
end
if isfield(Prob,'KNITRO')
   gnProb.KNITRO         = Prob.KNITRO;
end
gnProb.CGO.idea          = idea;
gnProb.CGO.epsRank       = epsRank;
gnProb.PriLevOpt         = DebugPriLev;
gnProb.CGO.backupSolver1 = backupSolver1;
gnProb.CGO.backupSolver2 = backupSolver2;


% Take advantage of any derivatives in local search

if isempty(dc)
   DerLvl_sn             = 1;
   DerLvl_gn             = 0;
   cDiff                 = 6;
else
   DerLvl_sn             = 3;
   DerLvl_gn             = 2;
   cDiff                 = 0;
end

snProb.GO.maxFunc1    = maxFunc1;
snProb.GO.maxFunc2    = maxFunc2;
snProb.GO.maxFunc3    = maxFunc3;
snProb.GO.DIRECT      = DIRECT;

switch lower(localSolver)
 case {'snopt','npsol','nlssol','minos'}
   gnProb.NumDiff        = 6;
   gnProb.SOL.optPar(39) = DerLvl_gn;
   gnProb.ConsDiff       = cDiff;
   snProb.SOL.optPar(39) = DerLvl_sn;
   snProb.ConsDiff       = cDiff;
   %snProb.optPar(10)     = 1E-8;
   snProb.optPar(12)     = 1E-8; % Minor optimality tolerance
 case {'glccluster'}
   if isempty(GOlocalSolver)
      snProb.GO.localSolver = localSolver;
   else
      snProb.GO.localSolver = GOlocalSolver;
   end
 case {'glcdirect','glbdirect'}
   snProb.NumDiff        = 0;
   snProb.ConsDiff       = 0;
 case {'glcfast','glbfast'}
   snProb.NumDiff        = 0;
   snProb.ConsDiff       = 0;
 otherwise
   gnProb.NumDiff        = 1;
   gnProb.ConsDiff       = cDiff > 0;
   snProb.ConsDiff       = cDiff > 0;
end
gnProb.GO.maxFunc1    = maxFunc1;
gnProb.GO.maxFunc2    = maxFunc2;
gnProb.GO.maxFunc3    = maxFunc3;
gnProb.GO.DIRECT      = DIRECT;
switch lower(globalSolver)
 case {'glccluster'}
   if isempty(GOlocalSolver)
      gnProb.GO.localSolver = localSolver;
   else
      gnProb.GO.localSolver = GOlocalSolver;
   end
 case {'glcdirect','glbdirect'}
   gnProb.NumDiff        = 0;
   gnProb.ConsDiff       = 0;
 case {'glcfast','glbfast'}
   gnProb.NumDiff        = 0;
   gnProb.ConsDiff       = 0;
 otherwise
end

Iter       = 0;
Update     = 1;
fnStar     = Inf;
alphaXXX   = 0;

z = Fpen;
% Set infeasible points as infinity before check on minimum
z(Fpen-F >= 1E-14)=Inf;

% Set integer infeasible points as infinity before check on minimum
if ~isempty(IntVars)
   XX=X(IntVars,:);
   if length(IntVars) == 1
      v = XX~=round(XX);
   else
      v = ~all(XX==round(XX));
   end
   z(v) = Inf;
end

[fMin,fIdx] = min(z);

if isinf(fMin)
   % No feasible point found
   Feasible = 0;
   % Take the best infeasible point
   [fMin,fIdx] = min(Fpen);
else
   Feasible = 1;
end

% Best point found in unit space
x_min = X(:,fIdx(1)); 

% Best point found in original space
if SCALE
   O_min = tomsol(9, x_L, x_min, x_D); 
else
   O_min = x_min;
end

FLOWER    = fIdx(1);
fMinIter  = 0;
fDiff     = Inf;
fDiffOld  = NaN;
NOUPDATE  = 0;
SAME1     = 0;
SAME2     = 0;
O_pre     = inf*ones(d,1);
alpha     = NaN;
% Find max Lipschitz constant in initial set
Dist      = tomsol(30,X,X);
LipMx     = 0;
for i = 1:size(X,2)-1
    LipMx = max(LipMx,max(abs(F(i)-F(i+1:end))./Dist(i+1:end,i)));
end

snProb.SIGN = 1;
gnProb.SIGN = 1;

% Set parameters used in gnSolve
gnProb.CGO.globalSolver = globalSolver;
gnProb.CGO.TargetMin    = TargetMin;
gnProb.PriLev           = PriSub;
gnProb.LOCAL            = LOCAL;
gnProb.epsX             = epsX;
if LOCAL
   gnProb.CGO.localSolver  = localSolver;
   gnProb.dLin             = dLin;
   gnProb.MaxIter          = MaxIter;
end
% Set parameters used in snSolve
gnProb.CGO.N            = N;         % Later copied to snProb
snProb.LOCAL            = LOCAL;
snProb.epsX             = epsX;
snProb.dLin             = dLin;
snProb.MaxIter          = MaxIter;
if isempty(GOlocalSolver)
   if strcmpi(globalSolver,localSolver)
      GOlocalSolver = GetSolver('con',1,0);
   else
      GOlocalSolver = localSolver; 
   end
end

% Becase new sn_f,sn_g tests on CGO.fnStar
gnProb.CGO.fnStar = NaN;
if isempty(IntVars)
   snProb = ProbCheck(snProb,localSolver,checkType('con'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('glc'));
else
   snProb = ProbCheck(snProb,localSolver,checkType('minlp'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('minlp'));
end

% *********************************************************
% *************** START MAIN ITERATION LOOP ***************
% *********************************************************

% n     = Number of points used in the interpolation
% nFunc = Total number of costly function evaluations
% nInit = Number of points used for startup, now or before
% Iter  = Number of fresh sampled points 
% modN  = Cycle step

modN       = -1-infStep;
addlocStep = 0; 
if isempty(x_opt)
   snOptf = []; 
   snOptg = [];
   snOptH = [];
   snOptE = [];
   dXO    = [];
else
   % Compute surface value/gradient/Hessian for optimal point, if given
   % global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];
   snProb.CGO.fnStar = NaN;
   snOptf = nlp_f(xOptS,snProb);
   snOptg = nlp_g(xOptS,snProb);
   snOptH = nlp_H(xOptS,snProb);
   snOptE = eig(snOptH);
   zz     = tomsol(30,xOptS,X);
   % Compute closest point to global optimum
   dXO    = min(zz);
   % Compute sum of all distances from global optimum to initial set X
   SumXO  = sum(zz);
   % Compute mean distance from global optimum to initial set X
   MeanXO = mean(zz);
end

Its  = [];
xInf = NaN*ones(d,1);
fInf = NaN;
if PriLev > 1 | IterPrint
   fprintf('Iter %3d n %3d ', 0, n);
      fprintf('nFunc %3d ', nFunc);
      tt = time([3 2 4 5]);
      if tt(4) < 10
         fprintf('%d/%d %d:0%1d ', tt);
      else
         fprintf('%d/%d %d:%2d ', tt);
      end
      fprintf('Cycle %2d ', modN);
      fprintf('               ');
      if ~isempty(fGoal) & ~isinf(fGoal)
         fprintf(' fGoal %8.5f', fGoal);
      end
      if Feasible
         fprintf(' fMinF %11.8f', fMin);
      else
         fprintf(' fMinI %11.8f', fMin);
      end
      maxF = max(F);
      fprintf(' max(F) %9.5g', maxF);
      fprintf(' median(F) %9.5g', median(F));
      fprintf(' range(F) %9.5g', maxF-fMin);

      %if TRANSFORM ~= 0
      %   fprintf('yMin %8.5f yNew %8.5f ', yMin, yNew);
      %end
      %fprintf('fE %8.5f ',fk);
      %fprintf('I%%');
      %fprintf(' %6.4f ',100*abs(EGOResult.f_k/yMin));
      fprintf('\n');
      xprint(O_min,'xMin:',' %12.8f',8)
      if ~isempty(x_opt)
         if PriLev > 2
            xprint(x_opt,'xOpt:',' %12.8f',8)
            fprintf('  SumXO (sum||xOpt-X||) %f ',SumXO);
            fprintf('MeanXO (mean||xOpt-X||) %f ',MeanXO);
            doO = min(tomsol(30,x_min,xOptS));
            fprintf(' doO (||xOpt-xMin|| %f ',doO);
            fprintf('\n');
            if PriLev > 3
               fprintf('  dXO (min||xOpt-X||) %f ',dXO);
               fprintf('  snOptf %8.5f ',snOptf);
               fprintf('snOptg: ');
               fprintf('%f ',snOptg);
               if length(snOptg) > 6, fprintf('\n'); end
               xprint(snOptE,'snEig: ');
               %disp(snOptH)
            end
         end
      else
          Its.dXO    = [];
          Its.snOptf = [];
          Its.snOptg = [];
          Its.snOpte = [];
      end
end

% -------------- Result saving -------------------

TIME0    = Prob.TIME0;
DIGIT    = 0.1;
TESTDIG  = [];
convflag = 0;
cpumax   = 0;
progress = 1;
FLOW     = nFunc;

snProb.ixDist = [];
gnProb.ixDist = [];
PARALLEL=0;
if PARALLEL
   xPARA  = [];
   fPARA  = [];
   F0PARA = F;
end

while nFunc < MaxFunc & convflag == 0
   
   if nFunc-FLOW > MaxCycle*(N+1), progress = 0; break; end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end

   % Set parameters in global structure CGO
   gnProb.CGO.X          = X;
   % gnProb.CGO.n        = n;
   % gnProb.CGO.d        = d;

   gnProb.CGO.rbfType    = rbfType;

   % Minimize s_n(y) using global optimization
   %  -  Solve  min[s_n(y)] with a local optimizer by starting from
   %     the interpolation point with the least function value

   Iter              = Iter + 1;

   % Set parameters used in snSolve
   snProb.snP              = Iter;
   snProb.x_0              = x_min;
   snProb.CGO              = gnProb.CGO;
   snProb.CGO.globalSolver = globalSolver;
   snProb.CGO.localSolver  = localSolver;
   snProb.CGO.fnStar       = NaN;
   % modN not yet updated, make preliminary update for snSolve
   snProb.CGO.modN         = mod(modN+1,N+1);
   snProb.PriLev           = PriSub;
   snResult                = snSolve(snProb);
   min_sn                  = snResult.f_k;
   if Feasible & min_sn-fMin > 1E-6*max(1,abs(fMin))
      fprintf('WARNING min_sn %15.12f > fMin %15.12f. ',min_sn,fMin);
      if ~isempty(backupSolver1)
         fprintf('Rerun with %s\n',backupSolver1);
         snProb.CGO.globalSolver = backupSolver1;
         snProb.CGO.localSolver  = backupSolver1;
         if strcmpi(backupSolver1,'oqnlp')
            snProb.OQNLP.options.RANDOMSEED = 12345;
         end
         snResult1               = snSolve(snProb);
         min_sn1                 = snResult1.f_k;
         snProb.CGO.globalSolver = globalSolver;
         snProb.CGO.localSolver  = localSolver;
         fprintf('%s new solution min_sn %15.12f\n',backupSolver1,min_sn1);
         if snResult1.ExitFlag == 0 & min_sn1 < min_sn
            min_sn   = min_sn1;
            snResult = snResult1;
         end
      end
      if min_sn-fMin > 1E-6*max(1,abs(fMin))
         fprintf('WARNING Still min_sn %15.12f > fMin %15.12f. ',min_sn,fMin);
         if ~isempty(backupSolver2)
            fprintf('Rerun with %s\n',backupSolver2);
            snProb.CGO.globalSolver = backupSolver2;
            snProb.CGO.localSolver  = backupSolver2;
            if strcmpi(backupSolver2,'oqnlp')
               snProb.OQNLP.options.RANDOMSEED = 54321;
            end
            snResult2               = snSolve(snProb);
            min_sn2                 = snResult2.f_k;
            snProb.CGO.globalSolver = globalSolver;
            snProb.CGO.localSolver  = localSolver;
            fprintf('%s new solution min_sn %15.12f\n',backupSolver2,min_sn2);
            if snResult2.ExitFlag == 0 & min_sn2 < min_sn
               min_sn   = min_sn2;
               snResult = snResult2;
            end
         end
      end
   end
   if isempty(snResult.x_k)
      min_sn_y = x_min;
   else
      min_sn_y = snResult.x_k(:,1);
   end

   % Transform to original space   
   if SCALE
      min_snc = tomsol(9, x_L, min_sn_y, x_D); 
      c0      = tomsol(9, x_L, x_min, x_D); 
   else
      min_snc = min_sn_y; 
      c0      = x_min;
   end


   % Choose a target value fnStar or compute alpha (implicitly gives fnStar)

   modN = modN + 1;
   if modN > N
      if FLOWER < nFunc - N
         %modN = -1
         modN = 0-infStep;
      else
         modN = 0-infStep;
      end
      % Important, modN must be less than -1 for locSteps
      if addlocStep > 0
         modN = -1-addlocStep; 
      end
   end
   % Important - needed to to handle locSteps if no infStep 
   if modN == -1 & infStep == 0
      modN = 0;
   end
   ucOK = 1;
   if modN == -1
      alpha  =  inf;
      fnStar = -inf;
   elseif modN < -1
      % Steps to add interior points now initiated, reset addlocStep
      addlocStep  = 0;
      alpha       =  inf;
      fnStar      = -inf;

   elseif idea == 1 
      if modN == 0
         nmax = n;
      else
         if DeltaRule
            % Current idea used:
            nmax = min(n,max(2,nmax-floor((nFunc-nInit)/N)));
         else
            % Always use ALL points
            nmax = n;
         end
      end
      
      if OldStrat
         F_sort = sort(F_m);
         max_F  = F_sort(nmax);
      else
         if REPLACE > 1 & nmax > n/2
            max_F = median(F);
         elseif REPLACE == 1 & nmax > n/2
            max_F = median(F);
         elseif REPLACE == 0
            F_sort = sort(F);
            max_F = F_sort(nmax);
         else
            F_sort = sort(F_m);
            max_F  = F_sort(nmax);
         end
      end
      
      fOld = fnStar;
   % HKHX Special - use any of this for update of fDiff ???
      %max_F = max(F);
      %if fMin >= 0
      %   fDiff = max(1,fMin);
      %elseif max_F <= 0
      %   fDiff = 0.5*(max_F-fMin);
      %else
      %   fDiff = max(1,(abs(fMin)));
      %end

      %if f_Low > -1E299
      %   fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
      %else
      %   fDiff = max_F-min_sn;
      %end
      %fDiff = max(fDiff,max_F-fMin);
      % Use fMin instead of min_sn, because surface way be wild
      if f_Low > -1E299
         fDiff = max(0,min(max_F-fMin,min_sn-f_Low));
      else
         fDiff = max_F-fMin;
      end
      if fDiff <= 0
         %HKH New
         fDiff = max(0,max(F)-fMin);
         if fDiff <= 0
            fDiff = max(1E-4,max(F)-min(F));
         end
      end

      if fStarRule == 1
         fnStar = min_sn - ( (N-modN)/N )^2*fDiff;
         % Check that the above is exactly the following:
         %fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
      elseif fStarRule == 2
         fnStar = min_sn - ( (N-modN)/N )*fDiff;
      elseif fStarRule == 3
         if modN == 0
            fnStar = min_sn;
         elseif abs(min_sn) < 1E-12
            fnStar = min_sn - ( (N-modN)/N );
         else
            fnStar = min_sn - ( (N-modN)/N )*abs(min_sn);
         end
      else
         %if f_Low > -1E299
         %   fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
         %else
         %   fDiff = max_F-min_sn;
         %end
         % Use fMin instead of min_sn, because surface way be wild
         if f_Low > -1E299
            fDiff = max(0,min(max_F-fMin,min_sn-f_Low));
         else
            fDiff = max_F-fMin;
         end
         if fDiff <= 0
            %HKH New
            fDiff = max(0,max(F)-fMin);
            if fDiff <= 0
               fDiff = max(1E-4,max(F)-min(F));
            end
         end
         switch modN
         case 0
           fnStar = min_sn - fDiff;
         case 1
           fnStar = min_sn - 0.5*fDiff;
         case 2
           fnStar = min_sn - 1E-1*fDiff;
         case 3
           fnStar = min_sn - 1E-2*fDiff;
         case 4
           fnStar = min_sn - 1E-4*fDiff;
         case 5
           fnStar = min_sn;
         end
      end
      % fprintf('Median=max_F %20.10f \n',max_F)

      if REPLACE == 1 & ...
         abs(fOld - fnStar) < 1E-4*max(abs(fOld),abs(fnStar)) & (fStarRule < 3)
         disp('ALARM. fnStar nearly equal to fnStar in last step');
         fprintf('fnStar %18.12f ',fnStar);
         fprintf('fMin %18.12f ',fMin);
         fprintf('fDiff %18.12f ',fDiff);
         fprintf('max_F %18.12f ',max_F);
         fprintf('min_sn %18.12f ',min_sn);
         fprintf('\n');
         nmax = nmax - 1;
         % Create fnStar without REPLACE
         F_sort = sort(F);
         max_F = F_sort(nmax);
         fnStar = min_sn - ( (N-modN)/N )^2*(max_F-min_sn);
         fprintf('New fnStar %18.12f ',fnStar);
         fprintf('\n');
         % Check that the above is exactly the following:
         % fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
      end
       
      if modN == N
         % Unconstrained cycle step. But min_sn must be sufficiently lower
   
         if fMin == 0
            maxF = max(F);
            if min_sn >= -eps_sn*min(1,maxF) 
               if maxF == 0
                  fnStar = -0.01;
               else
                  fnStar = -0.01*min(1,maxF);
               end
               ucOK = 0; 
            end
         elseif min_sn >= fMin-eps_sn*abs(fMin)
            fnStar = min_sn - 1E-2*abs(fMin);
            ucOK = 0;
         end
      end
      gnProb.CGO.fnStar = fnStar;
   elseif idea == 2
      max_F = max(F_m);

      % Generalized the fDiff strategy
      %if f_Low > -1E299
      %   fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
      %else
      %   fDiff = max_F-min_sn;
      %end
      %fDiff = max(fDiff,max_F-fMin);
      % Change, because surface way be wild
      if f_Low > -1E299
         fDiff = max(0,min(max_F-fMin,min_sn-f_Low));
      else
         fDiff = max_F-fMin;
      end
      if fDiff <= 0
         %HKH New
         fDiff = max(0,max(F)-fMin);
         if fDiff <= 0
            fDiff = max(1E-4,max(F)-min(F));
         end
      end
      % Try improving search if failure during a whole cycle
      %if fDiffOld == fDiff & FLOWER < nFunc - N
      %   fDiff = fDiff*0.1;
      %end

      alphaPrev = alpha;
      if OldAlpha
       switch modN
          case 0
            alpha = fDiff/2;
            if alpha == 0
               alpha = 1;
            end
          case 1
            alpha = fDiff/2;
            if alpha == 0
               alpha = 0.5;
            end
          case 2
            alpha = min(fDiff/2,1);
            if alpha == 0
               alpha = 0.3;
            end
          case 3
            if fMin == 0
               %if (fMin-min_sn) <= eps_sn, ucOK = 0; end
               if abs(fMin-min_sn) <= eps_sn, ucOK = 0; end
            %elseif (fMin-min_sn) <= eps_sn*abs(fMin)
            elseif abs(fMin-min_sn) <= eps_sn*abs(fMin)
               ucOK = 0;
            end
            if ucOK
               alpha = 0;
            else
               alpha = min(fDiff/2,0.5);
               if alpha == 0
                  alpha = 0.1;
               end
            end
       end
      else
       switch modN
          case 0
            if abs(Update) == 1
               alpha = fDiff/2;
            elseif alphaPrev >0
               alpha = alphaPrev/2;
            else
               alpha = fDiff/2;
            end
            if alpha <= 0
               alpha = 1;
               %alpha = alphaXXX;
               %alpha = 20;
            end
          case 1
            if fMinIter == Iter-1 
               %???alpha = fDiff/2;
               alpha = fDiff/3;
            elseif fDiffOld == fDiff
               alpha = fDiff/4;
            elseif abs(Update) == 1
               %alpha = fDiff/2;
               %alpha = fDiff/3;
               alpha = fDiff/4;
            elseif alphaPrev >0
               alpha = alphaPrev/2;
            else
               alpha =  fDiff/2;
            end
            if alpha <= 0
               alpha = 0.5;
               %alpha = 0.5*alphaXXX;
               %alpha = 10;
            end
          case 2
            if fMinIter == Iter-1 
               %alpha = min(fDiff/2,1);
               % testa
               alpha = alphaPrev/2;
            elseif fDiffOld == fDiff
               if fDiff/2 > 1
                  alpha = 1;
               else
                  alpha = fDiff/4;
               end
            elseif abs(Update) == 1
               %alpha = min(fDiff/2,1);
               alpha = min(fDiff/8,1);
            else
               alpha = alpha/2;
            end
            if alpha <= 0
               alpha = 0.3;
               %alpha = 0.3*alphaXXX;
               %alpha = 5;
            end
          case 3
            if fMin == 0
               %URKif (fMin-min_sn) <= eps_sn, ucOK = 0; end
               if abs(fMin-min_sn) <= eps_sn, ucOK = 0; end
            %URKelseif (fMin-min_sn) <= eps_sn*abs(fMin)
            elseif abs(fMin-min_sn) <= eps_sn*abs(fMin)
               ucOK = 0;
            end
            if ucOK
               alpha = 0;
            else
               if fMinIter == Iter-1 
                  %alpha = min(fDiff/2,1);
                  % testa
                  alpha = alphaPrev/2;
               elseif fDiffOld == fDiff
                   if fDiff/2 > 1
                      alpha = 0.1;
                   else
                      %alpha = min(fDiff/4,0.1);
                      alpha = min(fDiff/8,0.1);
                   end
               elseif abs(Update) == 1
                   %alpha = min(fDiff/4,0.1);
                   alpha = min(fDiff/8,0.1);
               else
                   alpha = alpha/2;
               end
            end
       end
      end
      if ~abs(Update) == 1 & alphaPrev == alpha
         if alpha <= 0
            alpha = 1;
         else
            alpha = 2*alpha;
         end
         if PriLev > -1
            fprintf('Alpha same as previous step, adjust to %10.3f\n',alpha);
         end
      end
      gnProb.CGO.alpha = alpha;
   end

   % No point in not accepting local solution if Pure Integer problem
   if length(IntVars) == d, ucOK = 1; end
   
   % ********** MINIMIZE gn(y) **********
   if modN < -1
      % AddSurfMin option - add one additional interior surface point
      xLoc  = xOL(:,1);
      fLoc  = fOL(1);
      if PriLev > 3
         fprintf('Use the interior point xLoc with surface f(x) %f\n',fLoc);
         disp(xLoc);
      end
      % Remove the current interior point from the list of interior points in (xOL,fOL)
      xOL = xOL(:,2:end);
      fOL = fOL(2:end);
      % Not set (xInf,fInf), use same as last iteration
      %xInf = xLoc;
      %fInf = Inf;
   end
   if (modN == N) & ucOK
      if SolveInf
         NLP_x=[]; NLP_f=[]; NARG = [];
         gnProb.CGO.modN   = -1;
         gnProb.CGO.fnStar = 0;
         gnProb.CGO.alpha  = 0;
         gnProb.PriLev     = PriSub; 
         gnProb.x_0        = min_sn_y;
         gnProb.X0         = snResult.multiMin.xOpt;
         gnResult          = gnSolve(gnProb);
         xInf              = gnResult.xLoc;
         fInf              = gnResult.fLoc;
         gnProb.CGO.modN   = modN;
         gnProb.CGO.fnStar = fnStar;
         gnProb.CGO.alpha  = alpha;
         NLP_x=[]; NLP_f=[]; NARG = [];
      end
      if PriLev > 3
         fprintf('Local search OK, no need switch to target value minimization of gn_f\n')
      end
      if AddSurfMin > 0
         % Pick up global multiMin solutions for analysis of option AddSurfMin
         IX   = snResult.multiMin.IX;
         xOpt = snResult.multiMin.xOpt;
         fOpt = snResult.multiMin.fOpt;
      end
      OLDuc = 0;
      if OLDuc
         xGlob = [];
         fGlob = [];
         xLoc  = min_sn_y;
         fLoc  = min_sn;
      else
         % Minimize s_n(y) using local optimization
         %  -  Solve  min[s_n(y)] with a local optimizer by starting from
         %     the interpolation point with the least function value
         % i.e. Find a local solution close to (x_min,fMin)
         % Set parameters used in snSolve
         snProb.CGO.globalSolver = GOlocalSolver;
         snProb.CGO.localSolver  = GOlocalSolver;
         snProb.CGO.fnStar       = NaN;
         snProb.x_0              = x_min;
         snProb.PriLev           = PriSub;
         snResult                = snSolve(snProb);
         fLoc                    = snResult.f_k;
         xLoc                    = snResult.x_k(:,1);
         %xprint(fLoc,'fLoc  ');
         %xprint(min_sn,'min_sn');
         %xprint(xLoc,'xLoc    ');
         %xprint(min_sn_y,'min_sn_y');
         Dist1 = tomsol(30,xLoc,x_min);
         Dist2 = tomsol(30,min_sn_y,x_min);
         if PriLev > 1 & Dist1 > Dist2
            fprintf('   USE global minimum min_sn_y instead!!! ')
            fprintf('||xLoc-xMin|| %f ', Dist1);
            fprintf('||min_sn_y-xMin|| %f ', Dist2);
            fprintf('\n');
            xLoc  = min_sn_y;
            fLoc  = min_sn;
         end
      end
      if AddSurfMin > 0
         ix      = [];
         LipSn   = [];
         LipfM   = [];
         LipQ    = [];
         for i = 1:length(IX)
             jj  = IX(i);
             xT  = xOpt(:,jj);
             d1 = norm(xT-xLoc);
             if PriLev > 3
                xprint(xT,'xTry:')
                fprintf('Index IX of xTry: %3d. ',jj);
                fprintf('||xTry-xLoc|| %f ',d1);
	     end
             if d1 > 0.01
                doTX = tomsol(30,xT,X);
                [doT,k]  = min(doTX);
                if PriLev > 3
                   fprintf('min_i ||xTry - X(:,i)|| %f ',doT);
	        end
                if doT > 0.01 % Accept
                   LipSn = abs(F(k)-fOpt(jj))/doT;
                   LipfM = abs(F(k)-fMin)/doT;
                   LipfMSn = LipfM/LipSn;
                   if PriLev > 3
                      fprintf('LipSn %8.2g ',LipSn);
                      fprintf('LipfM %8.2g ',LipfM);
                      fprintf('LipMx %8.2g ',LipMx);
                      fprintf('LipQ %8.2f ',LipfMSn);
	           end
                   % Check previously accepted points
                   if isempty(ix)
                      distIX = Inf;
                   else
                      distIX = min(tomsol(30,xT,xOpt(:,ix)));
                      if PriLev > 3
                         fprintf('distIX %8.2f ',distIX)
                      end
                   end
                   if LipfM <= LipMx & distIX > 0.01
                   % Accept interior point if Lipschitz constant < max Lipschitz
                   % and not too close to previously accepted interior points
                      ix  = [ix;jj];
                      LipQ(length(ix)) = LipfMSn;
                      if PriLev > 3
                         fprintf(' OK')
                      end
                   end
                end
             end
             if PriLev > 3
                fprintf('\n');
             end
         end
         if ~isempty(ix)
            % Check the accepted points
            if AddSurfMin < length(ix)
               % Must rank the minima, use the quote
               [vv,iy] = sort(LipQ);
               iz = ix(iy(1:AddSurfMin));
               xOL = xOpt(:,iz);
               fOL = fOpt(iz);
               if PriLev > 3
                  fprintf('Reduce the Accepted minima: ');
                  disp(xOL)
               end
            else
               xOL = xOpt(:,ix);
               fOL = fOpt(ix);
               if PriLev > 3
                  fprintf('Accept the minima: ');
                  xprinti(ix);
                  disp(xOL);
               end
            end
            addlocStep = length(fOL);
         end
      end
   end
   if modN >= -1 & (modN ~= N | (modN == N & ~ucOK)) 
      
      if modN == -1 | SolveInf
         gnProb.CGO.fnStar = 0;
         gnProb.CGO.alpha  = 0;
         NLP_x=[]; NLP_f=[]; NARG = [];
         gnProb.CGO.modN   = -1;
         gnProb.CGO.fnStar = 0;
         gnProb.CGO.alpha  = 0;
         gnProb.IX         = [];
         gnProb.xOpt       = [];
         gnProb.fOpt       = [];
         gnProb.PriLev     = PriSub; 
         gnProb.x_0        = min_sn_y;
         gnProb.X0         = snResult.multiMin.xOpt;
         gnResult          = gnSolve(gnProb);
         xLoc              = gnResult.xLoc;
         fLoc              = gnResult.fLoc;
         xInf              = xLoc;
         fInf              = fLoc;
         gnProb.CGO.modN   = modN;
         gnProb.CGO.fnStar = fnStar;
         gnProb.CGO.alpha  = alpha;
         NLP_x=[]; NLP_f=[]; NARG = [];
      end
      if modN ~= -1
         gnProb.CGO.modN = modN; 
         if fnStar > fMin & PriLev > 0
            fprintf('------ WARNING!!! fnStar > fMin)');
            fprintf('fnStar %18.12f ',fnStar);
            fprintf('fMin %18.12f ',fMin);
            fprintf('\n');
         end
         % Set iteration dependent parameters used in gnSolve, output from snSolve
         if TargetMin == 2
            % Save interior point solutions for possible use in gnSolve
            IX   = snResult.multiMin.IX;
            xOpt = snResult.multiMin.xOpt;
            fOpt = snResult.multiMin.fOpt;
         end
         gnProb.IX   = snResult.multiMin.IX;
         gnProb.xOpt = snResult.multiMin.xOpt;
         gnProb.fOpt = snResult.multiMin.fOpt;
         fnStar = gnProb.CGO.fnStar;

         if idea == 1 & modN == 0 
            onBMin = d+1;
            onBIx  = 1;
            AmpFac = Inf;
            if PriSub > 2 & PriLev > 1
               fprintf('\n--- NEW major cycle %d ---\n',modN);
	    end
            %vJ = [100 50 20 10 7 5 2 1 0.7 0.5 0.2 0.1 0.05 0.01 0.001 1E-4];
            vJ = [100 50 20 10 7 5 2 1 ];

%vJ  = [0 1E-4 2.5E-4 5E-4 7.5E-4 1E-3, 2.5E-3 5E-3 7.5E-3, 0.01*[1:13],0.15,0.20,0.25 0.30 0.40 0.50 0.75 1.00 1.5 2 3 1E2 -inf];
            for i=1:length(vJ)
               AmpTry            = vJ(i);
               fnStar            = min_sn - AmpTry*fDiff;
               gnProb.CGO.fnStar = fnStar;
               gnProb.PriLev     = PriSub; 
               if i == 1
                  gnProb.x_0        = min_sn_y;
                  gnProb.X0         = snResult.multiMin.xOpt;
               else
                  gnProb.x_0        = xTry;
                  gnProb.X0         = snResult.multiMin.xOpt;
               end
               gnTry             = gnSolve(gnProb);
               xTry              = gnTry.xLoc;
               fTry              = gnTry.fLoc;
               if isempty(IV)
                  onB = nOnBound(xTry,x_L,x_U,epsX);
               else
                  Reals = find(~IV);
                  onB = nOnBound(xTry(Reals),x_L(Reals),x_U(Reals),epsX);
               end
               if PriLev > 1
                  fprintf('Try %2d: ',i);
                  fprintf('modN %d ',modN);
                  fprintf('AmpFac %10.4f ',AmpTry);
                  fprintf('fnStar %12.5g ',fnStar);
                  fprintf('fDiff %12.5g ',fDiff);
                  fprintf('onB %d',onB);
               end
               if onB < onBMin
                  onBMin   = onB;
                  onBIx    = i;
                  AmpFac   = AmpTry;
                  gnResult = gnTry;
                  xLoc     = xTry;
                  fLoc     = fTry;
                  if PriLev > 1
                     fprintf(' Save!!!');
                  end
               end
               if PriLev > 1
                  fprintf('\n');
               end
               if onB == 0
                  break;
               end
            end
            vJ = [100 50 20 10 7 5 2 1 0.7 0.5 0.2 0.1 0.05 0.01 0.001 1E-4];
         elseif idea == 1 & modN > 0 & modN < N
            if PriSub > 2 & PriLev > 1
               fprintf('--- Next cycle %d: ',modN);
	    end
            onBMin = d+1;
            onBIx  = 1;
            AmpFac = Inf;
            vJ = [vJ, vJ(end)/2];
   
            %for i=max(5+2*modN,onBIx+1):length(vJ)
            %for i=5+2*modN:length(vJ)
            for i=5+2*modN:12+modN
               AmpTry            = vJ(i);
               fnStar            = min_sn - AmpTry*fDiff;
               gnProb.CGO.fnStar = fnStar;
               gnProb.PriLev     = PriSub; 
               if i == 1
                  gnProb.x_0        = min_sn_y;
                  gnProb.X0         = snResult.multiMin.xOpt;
               else
                  gnProb.x_0        = xTry;
                  gnProb.X0         = snResult.multiMin.xOpt;
               end
               gnTry             = gnSolve(gnProb);
               xTry              = gnTry.xLoc;
               fTry              = gnTry.fLoc;
               if isempty(IV)
                  onB = nOnBound(xTry,x_L,x_U,epsX);
               else
                  Reals = find(~IV);
                  onB = nOnBound(xTry(Reals),x_L(Reals),x_U(Reals),epsX);
               end
               if PriLev > 1
                  fprintf('Try %2d: ',i);
                  fprintf('modN %d ',modN);
                  fprintf('AmpFac %10.4f ',AmpTry);
                  fprintf('fnStar %12.5g ',fnStar);
                  fprintf('fDiff %12.5g ',fDiff);
                  fprintf('onB %d',onB);
               end
               if onB < onBMin
                  onBMin  = onB;
                  onBIx   = i;
                  AmpFac  = AmpTry;
                  gnResult = gnTry;
                  xLoc     = xTry;
                  fLoc     = fTry;
                  if PriLev > 1
                     fprintf(' Save!!!');
                  end
               end
               if PriLev > 1
                  fprintf('\n');
               end
               if onB == 0
                  break;
               end
            end
   
         else
            if PriSub > 2 & PriLev > 1
               fprintf('--- Other cycle %d: ',modN);
               fprintf('fnStar %12.5g',gnProb.CGO.fnStar);
               fprintf('\n');
	    end
            gnProb.PriLev  = PriSub; 
            gnProb.x_0     = min_sn_y;
            gnProb.X0      = snResult.multiMin.xOpt;
            gnResult       = gnSolve(gnProb);
            xLoc           = gnResult.xLoc;
            fLoc           = gnResult.fLoc;
         end
      end
   end
   
   % Best point found on surface is xLoc - set as xNew
   xNew = xLoc;
   ix = find(all(X==xNew*ones(1,size(X,2))));
   if ~isempty(ix)
      % New point is already sampled in a previous step
      RESCUE = 1;
      % Try to use minimum of the surface instead
      if ~(modN == N & ucOK)
         xNew  = min_sn_y;
         ix = find(all(X==xNew*ones(1,size(X,2))));
      end
   else
      RESCUE = 0;
   end
   % New point in original space
   if SCALE
      O_new   = tomsol(9, x_L, xNew, x_D); 
   else
      O_new   = xNew;
   end
   % Compute errors in the new point xNew (O_new)
   % global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];
   snNew = nlp_f(xNew,snProb);
   NLP_x=[]; NLP_f=[]; NARG = [];
   gnProb.CGO.modN   = -1;
   gnProb.CGO.fnStar = 0;
   gnProb.CGO.alpha  = 0;
   % my problem is using a -1/my objective function
   myNew = nlp_f(xNew,gnProb);
   if myNew ~= 0
      myNew = -1/myNew;
   else
      myNew = inf;
   end
   gnProb.CGO.modN   = modN;
   gnProb.CGO.fnStar = fnStar;
   gnProb.CGO.alpha  = alpha;
   
   
   % ********** PLOTTING **********
   if PLOT & d > 1 
      
      if 10 & ~(modN == N)
         switch lower(globalSolver)
         case 'glcdirect'
              glb=load('glcDirectSave.mat','C');
         case 'glbdirect'
              glb=load('glbDirectSave.mat','C');
         case 'glcfast'
              glb=load('glcFastSave.mat','C');
         case 'glbfast'
              glb=load('glbFastSave.mat','C');
         case 'glbsolve'
              glb=load('glbSave.mat','C');
         case 'glcsolve'
              glb=load('glcSave.mat','C');
         otherwise
              break;
         end
         % Transform to original space
         if SCALE
            CCC   = tomsol(9, x_L, glb.C, x_D); 
         else
            CCC   = glb.C;  
         end
         CCC = zeros(size(glb.C));
         
         plot(CCC(1,:),CCC(2,:),'.r');
         hold on
         plot(O(1,:),O(2,:),'*');
      else
         plot(O(1,:),O(2,:),'*');
         hold on
         TT = 0:0.1:1;
         YY(1,:) = c0(1)+TT*(O_new(1)-c0(1));
         YY(2,:) = c0(2)+TT*(O_new(2)-c0(2));
         plot(YY(1,:),YY(2,:),'-g')
      end   
      plot(O_new(1),O_new(2),'*g');
      x_opt = Prob.x_opt;
      if ~isempty(x_opt)
         if min(size(x_opt))==1
            x_opt = x_opt(:);
         end
         plot(x_opt(1,:),x_opt(2,:),'*k');
      end
      hold off
      grid on
      zoom on
   end % if plot
   if PLOT & d == 1 
         switch lower(globalSolver)
         case 'glcdirect'
              glb=load('glcDirectSave.mat','C');
         case 'glbdirect'
              glb=load('glbDirectSave.mat','C');
         case 'glcfast'
              glb=load('glcFastSave.mat','C');
         case 'glbfast'
              glb=load('glbFastSave.mat','C');
         case 'glbsolve'
              glb=load('glbSave.mat','C');
         case 'glcsolve'
              glb=load('glcSave.mat','C');
         otherwise
              break;
         end
         % Transform to original space
         if SCALE
            CCC   = tomsol(9, x_L, glb.C, x_D); 
         else
            CCC   = glb.C;  
         end
         CCC = zeros(size(glb.C));
      plot(CCC(1,:),zeros(size(CCC,2)),'.r');
      plot(O_new(1),fLoc,'*g');
      x_opt = Prob.x_opt;
      if ~isempty(x_opt)
         if min(size(x_opt))==1
            x_opt = x_opt(:);
         end
         plot(x_opt(1,:),x_opt(2,:),'*k');
      end
      hold off
      grid on
      zoom on
   end
   
   % *************** UPDATE ***************
   % Remove information for response surface problem
   % global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];

   % Distance between new point and best point found or closest point in X
   [onB, doX, doM] = statGN(gnProb.x_L,gnProb.x_U,xNew,x_min,X,epsX,Update);

   fRed = 0;
   if isempty(O_new) % Infeasibility problem
      fNew   = Inf;
      fPen   = Inf;
      SAME1  = 0;
      SAME2  = SAME2+1;
      Update = -5;
   else
      % CHECK if to reject or accept new point
      if all(O_new == O_min), SAME1 = SAME1 + 1; else SAME1 = 0; end
      if all(O_new == O_pre), SAME2 = SAME2 + 1; else SAME2 = 0; end
      if isempty(ix) 
         %if isempty(ix) & (n <= 2*d*(d+1) | (onB == 0 & modN == 1) | modN ~= 1)
         % if isempty(ix) & (n <= 2*d*(d+1) | ((onB == 0 & modN < N-1) | modN >= N-1))
         % if isempty(ix) & (n <= d*(d+1)/2 | ((onB == 0 & modN < N-2) | modN >= N-2))
         % if isempty(ix) & ((onB == 0 & modN < N) | modN == N)
         % if isempty(ix) & (n <= 2*d*(d+1)   | ((onB == 0 & modN < N) | modN == N))
         Update = 1;
         fNew  = nlp_f(O_new, Prob, varargin{:});
         nFunc = nFunc + 1; % Total number of sampled points
         fPen = fNew;
         if dLin > 0
            L = Prob.A*O_new;
            fPen = fPen+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
         end
         if dCon > 0 
            C = nlp_c(O_new, Prob, varargin{:});
            nCon = nCon + 1;
            fPen = fPen+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
         end
      else
         Update = -4;
         fNew   = NaN;
      end
   end
   time  = fix(clock);

   if PARALLEL
      if Update == 1
         xPARA = [xPARA,xNew];
         fPARA = [fPARA;fNew];
         if modN == N
            % Update RBF surface with all points in the cycle
            if REPLACE > 1
               for i=1:size(xPARA,2)
                   F0   = [F0PARA;fPara(1:i)];
                   FMIN = min(F0);
                   if FMIN <= 0
                      FMAX = 10^REPLACE;
                   else
                      FMAX = 10^(ceil(log10(FMIN))+REPLACE);
                   end
                   ix = find(F0>FMAX);
                   F_m0 = F0;
                   if ~isempty(ix)
                      F_m0(ix) = FMAX + log10(F0(ix)-FMAX+1);
                   end
                   Update = tomsol(24, xPARA(:,i), F_m0);    
               end
            elseif REPLACE == 1
               for i=1:size(xPARA,2)
                   fnewV=min(median([F0PARA;fPARA(1:i)]),[F0PARA;fPARA(1:i)]);
                   Update = tomsol(24, xPARA(:,i), fnewV);    
               end
            else
               for i=1:size(xPARA,2)
                   Update = tomsol(24, xPARA(:,i), fNew(i));    
               end
            end
            xPARA  = [];
            fPARA  = [];
            F0PARA = [F0PARA;fPARA];
         end
      end
      
   else
      if REPLACE > 1 & Update == 1
         F0   = [F;fNew];
         FMIN = min(F0);
         if FMIN <= 0
            FMAX = 10^REPLACE;
         else
            FMAX = 10^(ceil(log10(FMIN))+REPLACE);
         end
         ix = find(F0>FMAX);
         F_m0 = F0;
         if ~isempty(ix)
            fprintf('Restriced f values, indices:')
            xprinti(ix);
            F_m0(ix) = FMAX + log10(F0(ix)-FMAX+1);
         end
         Update = tomsol(24, xNew, F_m0);    
      elseif REPLACE == 1 & Update == 1
         F_m0=min(median([F;fNew(1)]),[F;fNew(1)]);
         Update = tomsol(24, xNew, F_m0);    
      elseif  Update == 1
         Update = tomsol(24, xNew, fNew);    
      end
   end
     
   if Update == 1 | Update == -1
      % All is OK, Update =-1, refactorization made the trick
      % New maximal Lipschitz constant
      LipMx = max(LipMx,max(abs(fNew-F)./tomsol(30,xNew,X)'));
      if REPLACE > 1
         F   = F0;
         F_m = F_m0;
      else
         F   = [F;fNew];
         if REPLACE == 1
            F_m = F_m0;
         else
            F_m = F; % No replacement
         end
      end
      Fpen= [Fpen;fPen];
      X   = [X xNew];
      O   = [O O_new];
      O_pre = O_new;
      n   = n + 1;
      
      % Is the function value in the new point less than fMin?
      % If feasible, compare fNew,fMin, otherwise fPen and fMin
      % Or now feasible, but previously not?
      fRed  = fMin - fPen; 
      if (Feasible & fNew < fMin & fPen==fNew) | ...
         (~Feasible & fPen < fMin) | (~Feasible & fPen == fNew)
         Feasible = fNew == fPen;
         fMin     = fPen;
         fMinIter = Iter;
         fIdx     = length(Fpen);
         x_min    = xNew;
         O_min    = O_new;
         FLOWER   = nFunc;
         FLOW     = nFunc;
         alphaXXX = alpha;
      end
      NOUPDATE = 0;
   elseif Update == -2
      % Update == -2 New point bad, even refactorization failed
      % Infeasible problem, ill-conditioning?

      Dist = tomsol(30,X(:,fIdx),X);
      Dist(fIdx) = Inf;

      control = -1;
      VALUE = 1;
      while control < 0
         tomsol(25) % Deallocates memory
         fprintf('Minimal distance to X set');
         [minDist ixDist] = min(Dist);
         fprintf(' %s',minDist);
         fprintf('\n');
         ix = ones(n,1);
         ix(ixDist)=VALUE;
         % Remove most infeasible point
         [maxPen ixF] = max(Fpen-F);
         if maxPen == 0
            % If all points feasible, remove point with largest f(x)
            [maxPen ixF] = max(Fpen);
         end
         ix(ixF)=VALUE;
         ix   = find(ix);
         F    = F(ix);
         Fpen = Fpen(ix);
         X    = X(:,ix);
         O    = O(:,ix);
         n    = size(X,2);
         Dist = Dist(ix);
         fIdx = find(isinf(Dist));
         %if VALUE == 1
         %   REPLACE=1-REPLACE
         %end
         if REPLACE > 1
            FMIN = min(F);
            if FMIN <= 0
               FMAX = 10^REPLACE;
            else
               FMAX = 10^(ceil(log10(FMIN))+REPLACE);
            end
            ix = find(F>FMAX);
            F_m = F;
            if ~isempty(ix)
               F_m(ix) = FMAX + log10(F(ix)-FMAX+1);
            end
         elseif REPLACE == 1
            F_m = min(median(F),F);
         else
            F_m = F;
         end
         % TOMSOL INIT, send F_m to Fortran, not F
         'make init again'
         control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
         VALUE = 0;
      end
      if control < 0
          fprintf('New initial interpolation failed');
          tomsol(25) % Deallocates memory
          Result.ExitFlag = 1;
          Result.ExitText = 'New Initial interpolation failed';
          Result.CGO.Its  = [];
          Result.DIGIT    = [];
          Result          = endSolve(Prob,Result);
          return     %Something is really wrong
      end
      NOUPDATE = 0;
   elseif Update == -3 | Update == -4 | Update == -5
      % Update == -3 Point too close to old point
      % Update == -4 Point identical to old point
      % Update == -5 No feasible point found on surface
      NOUPDATE = NOUPDATE+1;
      if PriLev > 1
         fprintf('!!!!! Update impossible\n')
      end
   end
   %if FLOWER < nFunc - N
   %   % Trying to improve a nonworking strategy, by changing replacement
   %   REPLACE = 1 - REPLACE;
   %end
   
   [onB_sn, doX_sn, doM_sn] = statGN(snProb.x_L,snProb.x_U,min_sn_y,...
       x_min,X,epsX,Update);
   if ~isempty(x_opt)
      doO_sn = min(tomsol(30,min_sn_y,xOptS));
   else
      doO_sn = [];
   end
   if ~isempty(x_opt)
      doO = min(tomsol(30,xNew,xOptS));
   else
      doO = [];
   end
   % Distance between new point and min on RBF surface
   SoO = min(tomsol(30,xNew,min_sn_y));
   if isempty(fGoal) | isinf(fGoal)
      RelErr = NaN;
   else
      if fGoal == 0
          RelErr = fMin-fGoal;
      else
          RelErr = (fMin-fGoal)/abs(fGoal);
      end
   end
   snErr = fNew-snNew;

   Its.Iter(Iter)=Iter;
   Its.n(Iter)=n;
   Its.nFunc(Iter)=nFunc;
   Its.modN(Iter)=modN;
   if idea == 1
      Its.fnStar(Iter)=fnStar;
   else
      Its.fnStar(Iter)=alpha;
   end
   Its.fDiff(Iter)=fDiff;
   Its.fMin(Iter)=fMin;
   Its.fRed(Iter)=fRed;
   Its.RelErr(Iter)=RelErr;
   Its.snErr(Iter)=snErr;
   Its.FLOWER(Iter)=FLOWER;
   Its.fMinIter(Iter)=fMinIter;
   Its.fNew(Iter)=fNew;
   Its.minSn(Iter)=min_sn;
   Its.onBminSn(Iter)=onB_sn;
   Its.distminSn2X(Iter)=doX_sn;
   Its.distminSn2xMin(Iter)=doM_sn;
   Its.onBxNew(Iter)=onB;
   Its.distxNew2X(Iter)=doX;
   Its.distxNew2xMin(Iter)=doM;
   Its.distxNew2minSn(Iter)=SoO;
   if ~isempty(x_opt)
      % global NLP_x NLP_f NARG
      NLP_x=[]; NLP_f=[]; NARG = [];
      snOptf = nlp_f(xOptS,snProb);
      snOptg = nlp_g(xOptS,snProb);
      snOptH = nlp_H(xOptS,snProb);
      snOptE = eig(snOptH);
      dXO    = min(tomsol(30,xOptS,X));

      Its.snOptf(Iter)   = snOptf;
      Its.snOptg(:,Iter) = snOptg;
      Its.snOptE(:,Iter) = snOptE;
      Its.dXO(Iter)      = dXO;
   end

   % New variables saved
   Its.LipMx(Iter)     = LipMx;
   Its.minSn_y(:,Iter) = min_sn_y;
   Its.xInf(:,Iter)    = xInf;
   Its.fInf(Iter)      = fInf;
   Its.snNew(Iter)     = snNew;

   if PriLev > 1 | IterPrint
      fprintf('Iter %3d n %3d ', Iter, n);
      fprintf('nFunc %3d ', nFunc);
      tt = time([3 2 4 5]);
      if tt(4) < 10
         fprintf('%d/%d %d:0%1d  ', tt);
      else
         fprintf('%d/%d %d:%2d ', tt);
      end
      fprintf('Cycle %2d', modN);
      if RESCUE
         fprintf('R');
      else
         fprintf(' ');
      end
      if idea == 1
         fprintf(' fnStar%7.3f',fnStar);
      else
         fprintf(' alpha %7.3f',alpha);
      end
      if ~isempty(fGoal) & ~isinf(fGoal)
         fprintf(' fGoal %8.5f', fGoal);
      end
      if Feasible
         fprintf(' fMinF %11.8f ', fMin);
      else
         fprintf(' fMinI %11.8f ', fMin);
      end
      if fRed > 0
         fprintf('at %3d/IT %d fNew %11.8f ', FLOWER, fMinIter,fNew);
      else
         fprintf('at %3d/It %d fNew %11.8f ', FLOWER, fMinIter,fNew);
      end
      if ~isnan(RelErr)
         fprintf(' RelErr %10.6f', RelErr);
      end
      fprintf(' snErr %10.6f', snErr);
      %fprintf(' min_sn %11.8f', min_sn);
      fprintf(' fLoc %11.8f', fLoc);
      %NEWHKHfprintf(' %11.8f', fGlob);
      fprintf('\n');
      xprint(O_new,'xNew:',' %12.8f',8)
      if PriLev > 2
         fprintf('  min_sn %f: ',min_sn);
         fprintf('[%d] ',onB_sn);
         fprintf('doX %f ',doX_sn);
         fprintf('doM %f ',doM_sn);
         if ~isempty(x_opt)
            fprintf('doO %f ',doO_sn);
         end
         fprintf('xNew: [%d] ',onB);
         fprintf('doX %f ',doX);
         fprintf('doM %f ',doM);
         if ~isempty(x_opt)
            fprintf('doO %f ',doO);
         end
         fprintf('SoO %f ',SoO);
         fprintf('\n');
      end
      if PriLev > 3
         fprintf('  snNew-min_sn %f ',snNew-min_sn);
         if idea == 1 
            fprintf('snNew-fnStar %f ',snNew-fnStar);
         end
         fprintf('snNew-fNew %f ',snNew-fNew);
         myNew = myNew/gnProb.SIGN;
         fprintf('myNew %f ',myNew);
         fprintf('fRed %f. ',fRed)
         if modN ~= N % | (Rescue == 0 & ~ucOK)
            %fprintf('hn %f ',-gnProb.SIGN/gnResult.f_k);
            fprintf('hn %f ',-1/gnResult.f_k);
            if idea == 1 
               %fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+1/gnResult.f_k);
               fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+gnProb.SIGN/gnResult.f_k);
            end
         end
         fprintf('\n');
      end
      if ~isempty(x_opt) & PriLev > 3
         fprintf('  dXO %f ',dXO);
         fprintf('snOptf %8.5f ',snOptf);
         fprintf('snOptg: ');
         fprintf('%f ',snOptg);
         if length(snOptg) > 6, fprintf('\n  '); end
         xprint(snOptE,'snEig: ');
         %disp(snOptH)
      end
   end
   
   % ********** CONVERGENCE TEST **********
   if n >= nMax, convflag = 7; end
   if convflag == 0
      convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
   end
   if convflag == 0
      if NOUPDATE > N
         convflag = 4;
      elseif SAME1 > N
         convflag = 5;
      elseif SAME2 > N
         convflag = 6;
      end
   end
   while RelErr <= DIGIT & DIGIT > 1E-6
      k = abs(log10(DIGIT));
      if (IterPrint | PriLev > 0) & convflag == 0
         fprintf('---------- MINIMUM REACHED WITH %d DIGITS!!\n',k);
      end

      TESTDIG.Iter(k)    = Iter;
      TESTDIG.FuncEv(k)  = nFunc;
      TESTDIG.CPUtime(k) = cputime-TIME0;
      DIGIT = DIGIT*0.1;
   end

   if PLOT
      pause(1)
   end
   if abs(Update) == 1
      fDiffOld = fDiff;
   end

   % -------------- Result saving -------------------
   %saveIter;
end

% *******************************************************
% *************** END MAIN ITERATION LOOP ***************
% *******************************************************

% SAVE RESULTS

fMinIdx = fIdx(1);
rngState = rand('state'); 

try
   save('cgoSave.mat','Name','O','F','X','F_m','nInit','Fpen','fMinIdx', ...
        'rngState');
catch
   warning('Failed to save warm start information to cgoSave.mat');
   disp('Check that the file is writable and that ego is executing in a writable folder');
   disp('Warm start information is returned in Result.CGO.WarmStartInfo');
end

% Create struct output with warmstart information
Result.CGO.WarmStartInfo = struct(...
   'Name',Name,...
   'O',O,...
   'F',F,...
   'X',X,...
   'F_m',F_m,...
   'nInit',nInit,...
   'Fpen',Fpen,...
   'fMinIdx',fMinIdx,...
   'rngState',rngState);

Result.CGO.snProb = snProb;
Result.CGO.gnProb = gnProb;
Result.CGO.fGoal  = fGoal;
Result.CGO.Its    = Its;

% All points i with F(i)=f_min
if SCALE
   Result.x_k      = tomsol(9,x_L,X(:,find(F==fMin)),x_D);    
else
   Result.x_k      = X(:,find(F==fMin));
end
if isempty(Result.x_k)
   if SCALE
      Result.x_k      = tomsol(9,x_L,X(:,find(Fpen==fMin)),x_D);    
   else
      Result.x_k      = X(:,find(Fpen==fMin));
   end
   fprintf('Warning: Optimal point is found to be infeasible: ');
   fprintf('constraints are violated\n')
end

Result.f_k      = fMin;     % Best function value
cMin = [];
if dCon > 0
   for i=1:size(Result.x_k,2)
       cMin = [cMin,nlp_c(Result.x_k(:,i), Prob, varargin{:})];
       nCon = nCon + 1;
   end
end
if dLin > 0 & size(Result.x_k,2)==1
   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at best x_k
end
Result.c_k      = cMin;     % Constraint value at best x_k
Result.Iter     = Iter;     % Number of iterations
Result.FuncEv   = nFunc;
Result.GradEv   = -1;
Result.ConJacEv = -1;
Result.ConstrEv = nCon;
Result.SolverAlgorithm = [Result.SolverAlgorithm ...
                          '. Global solver ' globalSolver  ...
                          '. Local solver ' localSolver];
Result.ExitFlag = 0;
Result.ExitText = ['Tried ' num2str(nFunc) ' f(x), using ' ...
                    num2str(n) ', startup ' num2str(nInit)];
if cpumax & convflag == 0
   Result.Inform   = 9;
   Result.ExitText = [Result.ExitText '. Max CPU reached. '];
elseif ~progress & convflag == 0
   Result.Inform   = 8;
   Result.ExitText = [Result.ExitText '. No progress for ' ...
                      num2str(nFunc-FLOWER) ' function evaluations'];
elseif convflag == 7
   Result.Inform   = convflag;
   Result.ExitText = [Result.ExitText, '. All feasible integers tried'];
else
   Result.Inform   = convflag;
end
Result.DIGIT       = TESTDIG;
Result             = endSolve(Prob,Result);

% -------------- Result saving -------------------
%saveResult;

tomsol(25) % Deallocates memory

function convflag = isClose(fGoal,f,fTol,nFunc,Iter,PriLev)

convflag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convflag = 1;
elseif fGoal == 0
   %if abs(f-fGoal) < fTol
   if abs(f) < fTol
      convflag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convflag = 3;
end

if convflag > 0 & PriLev > 0 
   if convflag == 1
      fprintf('\n\nFunction value %f is less than fGoal %f \n',f,fGoal);
   elseif convflag == 2
      fprintf('\n\nError in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal),fTol);
   elseif convflag == 3
      fprintf('\n\nRelative error in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('Number of function evaluations:  %d\n',nFunc);
   fprintf('Number of iterations:            %d\n',Iter);
end

% MODIFICATION LOG
%
% 071008  hkh  DYNRBF created from rbfSolve, implementing new algorithm ideas 
% 071010  hkh  Use gnProb.x_0 and gnProb.X0 as additional initial points
% 080617  frhe Made REPLACE>1 transformation monotonic
% 080617  frhe log10 used in REPLACE>1 in accordance with help
