% arbfMIP is based on the ARBF algorithm presented in
% 1. Kenneth Holmstrom: "An adaptive radial basis algorithm (ARBF) for
%    expensive black-box global optimization", Journal of Global Optimization,
%    2008.
% 2. Kenneth Holmstrom, Nils-Hassan Quttineh and Marcus M. Edvall, "An adaptive
%    radial basis algorithm {(ARBF)} for expensive black-box mixed-integer
%    constrained global optimization", Optimization and Engineering, 2008.
%
% arbfMIP solves problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finit
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Some or all x may be integer valued as specified by other input variables
%
% f(x) are assumed to be a costly function
% c(x) are assumed to be cheaply computed
%
% Any set {j} of costly constraints can be treated by adding penalty terms to 
% the objective function in the following way609G
%      p(x) = f(x) + SUM_j w_j * max (0,c^j(x)-c_U^j, c_L^j-c^j(x)),
% where weighting parameters w_j > 0 have been added.
% The user then returns p(x) instead of f(x) to the CGO solver.
%
% Calling syntax:
%
% function Result = arbfMIP(Prob, varargin)
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say RBFF
%   FUNCS.c   The routine to compute the nonlinear constraint, say RBFC
%             A call to tomFiles.m or glcAssign.m sets these fields.
%   x_L       Lower bounds for each element in x. Must be finite.
%   x_U       Upper bounds for each element in x. Must be finite.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   WarmStart If true, >0, arbfMIP reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%             If Prob.CGO.WarmStartInfo has been defined through a call to
%             WarmDefGLOBAL, this field is used instead of the cgoSave.mat file
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   user      User field used to send information to low-level functions
%   PriLevOpt Print Level
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%   PriLevSub Print Level in subproblem solvers, see help in snSolve, gnSolve
%   f_Low     Lower bound on the optimal function value. If defined, used to
%             restrict the target values into interval [f_Low,min(surface)]
% --------------------------------------------
% optParam    Structure in Prob, Prob.optParam
% ---------------------------------------------
%             Defines optimization parameters. Fields used:
%  MaxFunc    Maximal number of costly function evaluations, default 300 (<=5000)
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then MaxFunc = MaxFunc + nFunc
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1. See OUTPUT PRINTING below
%  fGoal      Goal for function value, if empty or Inf not used
%  eps_f      Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%             See the output field maxTri.
%  bTol       Linear constraint tolerance
%  cTol       Nonlinear constraint tolerance
%  MaxIter    Maximal number of iterations used by the local optimization
%             solver the response surface in each step. Default 1000,
%             except for pure IP problems, then max(GO.MaxFunc, MaxIter)
%
% ------------------
% Fields in Prob.CGO
% ------------------
%
% Percent     Type of strategy to get the initial sampled values:
%
%              Percent |  Exp. Design                   |  ExD
%             ---------|--------------------------------|-------
%                      |  CORNER STRATEGIES:            |
%                  900 |  All Corners                   |    1
%                  997 |  x_L + x_U + adjacent corners  |    2
%                  998 |  x_U + adjacent corners        |    3
%                  999 |  x_L + adjacent corners        |    4
%                      |                                |
%                      |  DETERMINISTIC STRATEGIES:     |
%                    0 |  User given initial points     |    5
%                99-94 |  Use glcDirect                 |    6
%                      |                                |
%                      |  LATIN BASED SAMPLING:         |
%                    1 |  Maximin LHD 1-norm            |    7
%                    2 |  Maximin LHD 2-norm            |    8
%                    3 |  Maximin LHD Inf-norm          |    9
%                    4 |  Minimal Audze-Eglais          |   10
%                    5 |  Minimax LHD (only 2 dim)      |   11
%                    6 |  Latin Hypercube               |   12
%                    7 |  Orthogonal Samling            |   13
%                      |                                |
%                      |  RANDOM STRATEGIES: (pp in %)  |
%                  1pp |  Circle surrounding            |   14
%                  2pp |  Ellipsoid surrounding         |   15
%                  3pp |  Rectangle surrounding         |   16
%
%             Negative values of Percent result in Constrained versions
%             of the Exp. Design methods 7-16. It means that all points
%             sampled are feasible with respect to all given constraints.
%
% nSample     Number of sample points to be used in initial experimental
%             design. nSample is used differently dependent on Percent:
%
%                     (n)Sample:
%               ExD |     < 0      |    = 0     |    > 0     |      []
%             ------|--------------|------------|------------|-------------
%                 1 |     2^d      |            |            |
%                 6 | abs(n) iters |    d+1     | max{d+1,n} | (d+1)(d+2)/2
%                 7 |   LATIN(k)   |            |            |
%              9-11 |     d+1      |            |            |
%             ------|--------------|------------|------------|-------------
%
%             where LATIN = [21 21 33 41 51 65 65] and k = abs(nSample).
%
%             Otherwise nSample as input does not matter.
%
%             DESCRIPTION OF THE EXPERIMENTAL DESIGNS:
%
%             ExD 1, All Corners. Initial points is the corner points of the
%             box given by Prob.x_L and Prob.x_U. Generates 2^d points, which
%             results in too many points when the dimension is high.
%
%             ExD 2, Lower and Upper Corner point + adjacent points.
%             Initial points are 2*d + 2 corners: the lower left corner x_L
%             and its d adjacent corners x_L+(x_U(i)-x_L(i))*e_i, i=1,...,d
%             and the upper right corner x_U and its d adjacent corners
%             x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 3. Initial points are the upper right corner x_U and its
%             d adjacent corners x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 4. Initial points are the lower left corner x_L and its
%             d adjacent corners x_L + (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 5. User given initial points, given as a matrix in CGO.X.
%             Each column is one sampled point. If d = length(Prob.x_L),
%             then size(X,1) = d, size(X,2) >= d+1. CGO.F should be defined
%             as empty, or contain a vector of corresponding f(x) values.
%             Any CGO.F value set as NaN will be computed by solver routine.
%
%             ExD 6. Use determinstic global optimization methods to find the
%             initial design. Current methods available (all DIRECT methods).
%
%             Percent:   99 = glcDirect,    97 = glcSolve,    95 = glcFast
%                        98 = glbDirect,    96 = glbSolve,    94 = glbFast
%
%             ExD 7-11. Optimal Latin Hypercube Designs (LHD) with respect to
%             different norms. The following norms and designs are available:
%
%             Percent:    1 = Maximin 1-Norm,         4 = Audze-Eglais Norm
%                         2 = Maximin 2-Norm,         5 = Minimax 2-Norm
%                         3 = Maximin Inf-Norm,
%
%             All designs taken from:  http://www.spacefillingdesigns.nl/
%
%             Constrained versions will try bigger and bigger designs up to 
%             M = max{10*d,nTrial} different designs, stopping when it has
%             found nSample feasible points.
%
%             ExD 12. Latin hypercube space-filling design. For nSample < 0,
%             k = abs(nSample) should in principle be the problem dimension.
%             The number of points sampled is:
%                   k      :  2  3  4  5  6  >6
%                   Points :  21 33 41 51 65 65
%             The call made is: X = daceInit(abs(nSample),Prob.x_L,Prob.x_U);
%
%             Set nSample = []  to get (d+1)*(d+2)/2 sampled points:
%                   d      : 1  2  3  4  5  6  7  8  9 10
%                   Points : 3  6 10 15 21 28 36 45 55 66
%             This is a more efficient number of points to use.
%
%             If CGO.X is nonempty, these points are verified as in ExD 5,
%             and treated as already sampled points. Then nSample additional
%             points are sampled, restricted to be close to the given points.
%
%             Constrained version of Latin hypercube only keep points that
%             fulfill the linear and nonlinear constraints. The algorithm
%             will try up to M = max{10*d,nTrial} points, stopping when it
%             has found nSample feasible points (d+1 points if nSample < 0).
%
%             ExD 13. Orthogonal Sampling, LH with subspace density demands.
%
%             ExD 14,15 and 16. Random strategies, the abs(Percent) value gives
%             the percentage size of an ellipsoid, circle or rectangle around
%             the so far sampled points that new points are not allowed in.
%             Range 1%-50%. Recommended values 10% - 20%.
%
%             If CGO.X is nonempty, these points are verified as in ExD 5,
%             and treated as already sampled points. Then nSample additional
%             points are sampled, restricted to be close to the given points.
%
% X,F,CX:     The fields X,F,CX are used to define user given points
%             ExD = 5 (Percent = 0) needs this information
%             If ExD == 6-12,14-16 these points are included into the design.
%
% X           A matrix of initial x values. One column for every x value. 
%             If ExD == 5, size(X,2) >= dim(x)+1 needed.
% F           A vector of initial f(x) values. If any element is set
%             to NaN it will be computed.
% CX          Optionally a matrix of nonlinear constraint c(x) values.
%             If nonempty, then size(CX,2) == size(X,2). If any element
%             is set as NaN, the vector c(x) = CX(:,i) will be recomputed.
%
% RandState    If >=0,  rand('state',RandState) is set to initialize the
%                       pseudo-random generator
%              If < 0   rand('state',sum(100*clock)) is set to give a new
%              or = []  set of random values each run.
%              If isnan(RandState), the random state is not initialized.
%
%              Default RandState = 0
%              RandState will influence if a stochastic initial experimental
%              design is applied, see input Percent and nSample.
%              RandState will also influence if using the multiMin solver,
%              but the random state seed is not reset in multiMin.
%              The state of the random generator is saved in rngState, and the
%              random generator is reinitialized with this state if WarmStart
%
% AddMP        If = 1, add the midpoint as extra point in the corner strategies.
%              Default AddMP=1 if any corner strategy.
%
%
% nTrial       For CLH, the method generates M = max{10*d,nTrial} trial points,
%              and evaluate them until nSample feasible points are found.
%
%              In the random designs, nTrial is the maximum number of trial
%              points randomly generated for each new point to sample.
%
% CLHMethod    Different search strategies for finding feasible LH points.
%              First of all, the least infeasible point is added. Then the
%              linear feasible points are considered. If more points are
%              needed still, the nonlinear infeasible points are added.
%
%              1 - Take the sampled infeasible points in order.
%              2 - Take a random sample of the infeasible points.
%              3 - Use points with lowest cErr.
%
% SCALE        0 - Original search space
%              1 - Transform search space to unit cube (default).
%
% REPLACE      0-No replacement (Default)
%               1 - Large function values are replaced by the median
%              >1 - Large values Z are replaced by new values
%              Replacement: Z:= FMAX + log10(Z-FMAX), where
%              FMAX = 10^REPLACE, if min(F) < 0
%              FMAX = 10^(ceil(log10(min(F)))+REPLACE), if min(F) >= 0
%              New replacement in every iteration, because min(F) may change
%
% SMOOTH       =1 The problem is smooth enough for local search using
%              numerical gradient estimation methods
%              =0 The problem is nonsmooth or noisy, and local search methods
%              using numerical gradient estimation are likely to produce
%              garbage search directions.
%              Default 1.
%
% globalSolver Global optimization solver used for subproblem optimization.
%              Default glcCluster (SMOOTH=1) or glcDirect (SMOOTH=0)
%              If the globalSolver is glcCluster, the fields
%              Prob.GO.maxFunc1, Prob.GO.maxFunc2, Prob.GO.maxFunc3,
%              Prob.GO.localSolver, Prob.GO.DIRECT and other fields
%              set in Prob.GO are used.
%              See the help for these parameters in glcCluster.
% localSolver  Local optimization solver used for subproblem optimization.
%              If not defined, the TOMLAB default constrained NLP solver is used
%
%              --- Special ARBF algorithm parameters ---
%
% rbfType      Type of radial basis function.
%              1 - Thin Plate Spline, 2 - Cubic (Default).
%              3 - Multiquadric,      4 - Inverse multiquadric
%              5 - Gaussian,          6 - Linear
%
% idea         Global search type, always idea = 1, i.e. use fnStar values
%
% N            Cycle length in idea 1 (default N=5 for fStarRule 1 and 2,
%              default N=1 for fStarRule 3) or idea 2 (always N=3)
%
% infStep      If =1, add search step with target value -inf first in cycle
%              Default 0.
% TargetMin    Which minimum of several to pick in target value problem
%              =0 Use global minimum
%              =1 Use best interior local minima, if none use global minimum
%              =2 Use best interior local minima, if none use RBF interior min
%              =3 Use best minimum with lowest number of coefficients on bounds
% fStarRule    Global-Local search strategy in idea 1. N = cycle length
%              Define min_sn as global minimum on surface. fStar Target value
%              1: fStar = min_sn - ((N-(n-nInit))/N)^2*Delta_n (Default)
%              2: fStar = min_sn -  (N-(n-nInit))/N   *Delta_n
%              Strategy 1 and 2 depends on Delta_n estimate (see DeltaRule)
%              If infStep true, addition of -inf-step first in cycle
%              3: fStar = -inf-step, min_sn-k*0.1*|min_sn| k=N,...,0
%
%              Strategy names in Gutmanns thesis: III, II, I
% DeltaRule    1 = Skip large f(x) when computing f(x) interval Delta
%              0 = Use all points. Default 1.
%
% eps_sn       Relative tolerance used to test if minimum of surface, min_sn, is
%              sufficiently lower than the best point found (fMin).  Default = 1E-7
% MaxCycle     Max number of evaluations without progress before stopping, default 10
%              Used in 3 convergence tests, see output parameter Inform, values 4,5,6.
% usecgolib    0 - do not use CGOLIB
%              1 - use CGOLIB
%
% ---------------------------------------------------------
% Fields in Prob.GO (Default values are set for all fields)
% ---------------------------------------------------------
%
% MaxFunc      Maximal number of function evaluations in each global search
% MaxIter      Maximal number of iterations in each global search
%    The following 5 fields are used by glcCluster, see help glcCluster
% DIRECT       DIRECT solver used in glcCluster: glcSolve or glcDirect(Default)
% maxFunc1     Maximum function evaluations 1st call to DIRECT routine
% maxFunc2     Maximum function evaluations 2nd call to DIRECT routine
% maxFunc3     Maximum sum of function evaluations in repeated 1st calls
%              to DIRECT routine trying to get feasible
% localSolver  The local solver used by glcCluster. If not defined, then
%              Prob.CGO.localSolver is used
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
% --------------------------------------------------------------------
% varargin    Additional parameters to arbfMIP are sent to the costly f(x)
% --------------------------------------------------------------------
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
%  ExitFlag Always 0
%  Inform   0 = Normal termination
%           1 = Function value f(x) is less than fGoal
%           2 = Error in function value f(x), abs(f-fGoal) <= fTol, fGoal=0
%           3 = Relative Error in function value f(x) is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           4 = No new point sampled for MaxCycle iteration steps
%           5 = All sample points same as the best point for MaxCycle last iter
%           6 = All sample points same as previous point for MaxCycle last iter
%           7 = All feasible integers tried
%           9 = Max CPU Time reached
%
% To make a warm start possible, arbfMIP saves the following information in
% the file cgoSave.mat:
% Name      Name of the problem
% O         Matrix with sampled points (in original space)
% X         Matrix with sampled points (in unit space if SCALE == 1)
% F         Vector with function values (penalty added for costly Cc(x))
% F_m       Vector with function values (replaced)
% F00       Vector of pure function values, before penalties
% Cc        Matrix with costly constraint values, Cc(x)
% nInit     Number of initial points >= d+1 (2^d if center points)
% Fpen      Vector with function values + additional penalty if infeasible
%           using the linear constraints and noncostly nonlinear c(x) 
% fMinIdx   Index of the best point found
% rngState  Current state of the random number generator used
%
% If the cgoSave.mat file fails to open for writing, the information is also
% available in Result.CGO.WarmStartInfo. Through a call to WarmDefGLOBAL,
% the Prob structure can be setup for warm start. In this case, th solver will 
% not load the data from cgoSave.mat.
%
% ---------------------------------------------
% OUTPUT PRINTING (IterPrint == 1 or PriLev > 0)
% ---------------------------------------------
% --- Row 1
% Iter      Number of iterations
% n         Number of trial x, n-Iter is number of points in initial design
% nFunc     Number of costly f(x) computed, nFunc <= n, n-nFunc = rejected pnts
% --->>     Time stamp (date and exact time of this printout)
% Cycle     Cycle steps: global=0 local=1 . infStep is marked -1,
%           Negative values - other type of steps
% fnStar    Target value fn_star
% fGoal     Goal value (if set)
% fMin      Best f(x) found so far. E.g. at 27/It 12 means n=27, Iter=12
%           fMinI means the best f(x) is infeasible
%           fMinF means the best f(x) is feasible (also integer feasible)
%           IT implies reduction in last step, It no reduction last step
%
% -------------------------------------------
% Additional information in iteration 1,2,...
% -------------------------------------------
% fNew      f(xNew), the costly function value for point tried in current Iter
% RelErr    Relative distance to known or assumed global minimum (Prob.x_opt)
% ---------------------------------------------------
% Additional information in iteration 0, on next 3 lines
% ---------------------------------------------------
% max(F)    maximum of all f(x) in the initial set of points X
% med(F)    median of all f(x) in the initial set of points X
% rng(F)    maxF-fMin, the range of f(x) values in the initial set X
% pDist     The size of the simply bounded region, ||x_U-x_L||_2
% LipU      Maximum Lipschitz constant for initial set X
% LipL      Minimum Lipschitz constant for initial set X
%
% xMin      Best point in initial set X
%
% xOpt      User-given global optimum Prob.x_opt (if defined)
%
% -----------
% In iteration 0 (if global optimum known and given in Prob.x_opt):
% -----------
% dXO       Minimal distance from global optimum to closest point of all
%           sampled points X in experimental design
% SumXO     Sum of distances from global optimum to all
%           sampled points X in experimental design
% doO       Distance from xBest with best f(x) in experimental design
%           to global optimum
% -----------
% --- Row 2
% -----------
% xNew      Point tried in the current iteration, scaled back if SCALE==1
% -----------
% --- Row 3 (if PriLev > 2) (All distances are in SCALED space [0,1]^d )
% -----------
% snErr     surfErr=Costly f(x) value - Surface value at x (Actual-Predicted)
% fLoc      Optimal subproblem solution
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
% doS       Distance from min_sn_y to new point xNew
% -----------
% --- Row 4 (if PriLev > 3) ???
% -----------
% snNew-min_sn
% snNew-fnStar
% snNew-fNew
% myNew
% fRed
% hn
% hnErr
% -----------
% --- Row 4 (if PriLev > 2)
% -----------
% GlobRed
% LocRed
% fRed
% LipUpp    Maximum Lipschitz constant for initial set X
% LipLow    Minimum Lipschitz constant for initial set X
%
% ======
% USAGE:
% ======
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
%    Prob.optParam.MaxFunc = 400;   % Change max number of function evaluations
%    Prob.optParam.MaxIter = 2000;  % Change local iteration maximum to 2000
%    Prob.GO.MaxFunc       = 20000; % 20000 global function values in subproblem
%
% Direct solver call:
%    Result = arbfMIP(Prob);
%    PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('arbfMIP',Prob,2);
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
% To make a restart, just set the restart flag, and call arbfMIP once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('arbfMIP',Prob,2);
%
% Another option for warm start is to put the restart information from the
% Result structure to the Prob structure and then call a CGO solver again:
%    Prob = WarmDefGLOBAL('arbfMIP',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('arbfMIP',Prob,2);
%
% It is also possible to run another CGO solver in the warm start:
%    Prob = WarmDefGLOBAL('arbfMIP',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('ego',Prob,2);
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 12, 2000. Last modified Oct 20, 2009.

% HKH unsolved issues
% AddSurfMin - use or not?
% varargin - possible to use?
% documentation of all output
% epsX set consistently
% Use of vjFac, increase - decrease, now limit to [1E-3,1E3]

function Result = arbfMIP(Prob, varargin)

if nargin < 1
   error('arbfMIP needs one parameter, the structure Prob');
end
global NLP_x NLP_f NARG

DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt
DebugARBF   = 0;  % If 1, use true f(x) values and compute f(x)
DP          = DebugARBF; % Short form of DebugARBF, used in assignments
DEBUG       = 0;  % DEBUG flag for tomsol
FASTCODE    = 1;  % Flag to avoid computing GO for all fnStar

time                   = fix(clock);
solvType               = checkType('glc');
Prob                   = ProbCheck(Prob,'arbfMIP',solvType);
Prob                   = iniSolve(Prob,solvType,0,0);

Result                 = ResultDef(Prob);
Result.Solver          = 'arbfMIP';
Result.SolverAlgorithm = 'Radial Basis Function Interpolation';

% Pick up bounds from Prob structure and check if OK
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds

if isempty(x_L) | isempty(x_U)
   disp('arbfMIP requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'arbfMIP requires both lower and upper variable bounds';
   Result.DIGIT    = [];
   Result          = endSolve(Prob,Result);
   return;
end
if ~(all(isfinite(x_L)) & all(isfinite(Prob.x_U)))
   disp('arbfMIP only solves box bounded problems, where both');
   disp('lower bounds Prob.x_L and upper bounds Prob.x_U are finite');
   Result.ExitFlag = 1;
   Result.ExitText = 'Problem not box bounded, variable bounds on x not finite';
   Result.DIGIT    = [];
   Result          = endSolve(Prob,Result);
   return;
end
% -----------------------------------------------------
% Get standard CGO variables from Prob.CGO structure
% and get variables from Prob.GO and Prob.MIP structures
% -----------------------------------------------------
[SCALE,REPLACE,RandState,Percent,nSample,AddMP,nTrial,CLHMethod,...
   SMOOTH,globalSolver,localSolver,LOCAL,GOMaxFunc,GOMaxIter,...
   maxFunc1,maxFunc2,maxFunc3,GOlocalSolver,DIRECT,IntVars,Reals,...
   nMax,d,x_L,x_U,x_D] = getCGOProbVars(Prob, x_L, x_U);
% -------------------------------------------------------------
% Set possibly redefined RandState, SCALE and IntVars into Prob
% -------------------------------------------------------------
Prob.CGO.RandState = RandState;
Prob.MIP.IntVars   = IntVars;
Prob.CGO.SCALE     = SCALE;

% HKH TEMP variables in new beta version
Prob.CGO.PEN     = 10;
% PEN=Prob.CGO.PEN; % PEN NOT USED NOW
SAVE = 1;

% ----------------------------------------------------------
% Set possibly redefined GOlocalSolver and DIRECT into Prob.GO
% ----------------------------------------------------------
Prob.GO.localSolver  = GOlocalSolver;
Prob.GO.DIRECT       = DIRECT;

% -------------------------------------------------
% Get special RBF variables from Prob.CGO structure
% -------------------------------------------------
[idea, rbfType, N, infStep, fStarRule, ...
   DeltaRule, MaxCycle, eps_sn, TargetMin, AddSurfMin, ...
   usecgolib] = getARBFProbVars(Prob, x_L, x_U);

% -------------------------------------------------------------
% Get variables from Prob structure
% Redefine GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3
% and set back in Prob.GO
% -------------------------------------------------------------
[Prob, MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, ...
   fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX, ...
   x_LL,x_UU, x_DD, pDist, dLin, dCon, ...
   Percent, nSample, nTrial, AddMP, REPLACE, ...
   backupSolver1, backupSolver2, xOptS ...
   ] = getProbVars(Prob, x_L,x_U,x_D, SCALE,globalSolver,IntVars,...
   Percent, nSample, nTrial, AddMP, REPLACE, ...
   GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3,'arbfMIP');


%  STEP 1, Initialization
%

% ------------------------------------------------------------
% Initial design from warm start or expDesign, set random seed
% ------------------------------------------------------------
if Prob.WarmStart
   % Restart with values from previous run
   [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, MaxFunc] = ...
      CGOWarmStart(Prob, MaxFunc, PriLev > 0 | IterPrint);
   if nFunc == 0
      Prob.WarmStart = 0;
   else
      ExDText = ['. Warm Start with ',num2str(nFunc),' f(x)'];
   end
   nCon = 0; % Count number of calls to constraint routine, nlp_c
end
if ~Prob.WarmStart
   [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, nCon, nSample, ...
      ExDText] = CGOInit( Prob, RandState, nSample, Percent, AddMP, nTrial,...
      CLHMethod, SCALE, REPLACE, PriLev > 0 | IterPrint, varargin);
end


% ----------------------------------------------------
% Save input parameters in output structure Result.CGO
% ----------------------------------------------------

switch rbfType
case 1
   Result.SolverAlgorithm = ['Thin-plate Spline ', Result.SolverAlgorithm];
case 2
   Result.SolverAlgorithm = ['Cubic ', Result.SolverAlgorithm];
case 3
   Result.SolverAlgorithm = ['Multiquadric ', Result.SolverAlgorithm];
case 4
   Result.SolverAlgorithm = ['Inverse Multiquadric ', Result.SolverAlgorithm];
case 5
   Result.SolverAlgorithm = ['Gaussian ', Result.SolverAlgorithm];
case 6
   Result.SolverAlgorithm = ['Linear ', Result.SolverAlgorithm];
end

Result.CGO = struct(...
   'SCALE',SCALE, 'REPLACE',REPLACE, 'RandState',RandState,...
   'Percent',Percent,'nSample',nSample,'AddMP',AddMP, ...
   'nTrial',nTrial, 'CLHMethod',CLHMethod, ...
   'fTol',fTol, 'f_Low',f_Low,'SMOOTH', SMOOTH, ...
   'idea',idea, 'rbfType',rbfType, 'N',N, 'infStep',infStep, ...
   'fStarRule',fStarRule, 'DeltaRule',DeltaRule, ...
   'MaxCycle',MaxCycle, 'eps_sn', eps_sn, ...
   'TargetMin',TargetMin, 'AddSurfMin',AddSurfMin, ...
   'globalSolver',globalSolver, 'localSolver',localSolver,...
   'LOCAL',LOCAL,...
   'GOMaxFunc',Prob.GO.MaxFunc,'GOMaxIter',Prob.GO.MaxIter,...
   'GOmaxFunc1',Prob.GO.maxFunc1,'GOmaxFunc2',Prob.GO.maxFunc2,...
   'GOmaxFunc3',Prob.GO.maxFunc3,...
   'backupSolver1',backupSolver1,'backupSolver2',backupSolver2);

Prob.CGOLIB = struct('usecgolib', usecgolib);

% ------------------------------------------------------------------------
% Initial RBF interpolation surface
% ------------------------------------------------------------------------
% TOMSOL INIT, send F_m to Fortran, not F
clear('tomsol');
% HKH Special safe guard now - to be investigated
MaxFunc=min(2000,MaxFunc);

if(usecgolib == 1)
  % Create the observation set
  Prob.CGOLIB.obs = cgolib(100, size(X,1), MaxFunc, 4);
  % Create the data column
  Prob.CGOLIB.datacol = cgolib(103, Prob.CGOLIB.obs, 1);
  % Add row data (original data)
  cgolib(102, Prob.CGOLIB.obs, X, Prob.CGOLIB.datacol, F);
  % If we are using the REPLACE transform, then create
  % a new column with this transform and base the RBF
  % surface on that column.
  
  % Transformation column doesn't exist yet.
  if(REPLACE == 1)
    Prob.CGOLIB.repcol = cgolib(111, Prob.CGOLIB.obs, 1, Prob.CGOLIB.datacol);
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.repcol;
    cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
  elseif(REPLACE == 0)
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.datacol;
  elseif(REPLACE > 1)
    Prob.CGOLIB.repcol = cgolib(111, Prob.CGOLIB.obs, 6, ...
				Prob.CGOLIB.datacol, ...
				{abs(REPLACE)});
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.repcol;
    cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
  else
    error('Invalid replace');
  end
  % Create RBF surface
  Prob.CGOLIB.rbf = cgolib(300, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
  % Choose RBF type
  % How?
  % Now interpolate
  cgolib(302, Prob.CGOLIB.rbf);
  % NEED CHECKS FOR THE INTERPOLATION!
else
  control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
  if control < 0
    fprintf('Initial interpolation failed');
    tomsol(25); % Deallocates memory
    Result.ExitFlag = 1;
    Result.ExitText = 'Initial interpolation failed';
    Result.DIGIT    = [];
    Result.CGO.Its  = [];
    Result          = endSolve(Prob,Result);
    return
  end
end

% -----------------------------------------
% Generate Prob structures for sub problems
% -----------------------------------------

%gProb = CGOGlobalProb(Prob, goName, glob_f, glob_g, glob_H,...
%   glob_c, glob_dc, glob_d2c, x_L, x_U, x_D, x_LL, x_UU, SCALE,...
%   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
%   backupSolver1, backupSolver2, PriSub)

snProb = CGOGlobalProb(Prob, 'RBFsn', 'sn_f', 'sn_g', [],...
   'rbf_c', 'rbf_dc', 'rbf_d2c', x_L, x_U, x_D, x_LL, x_UU, SCALE,...
   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
   backupSolver1, backupSolver2, PriSub);

snProb.CGOLIB = Prob.CGOLIB;
snProb.GO.ProbL.CGOLIB = Prob.CGOLIB;

% --------------------------------------------------
% Special input in snProb
% --------------------------------------------------
if DebugPriLev > 0
   snProb.optParam.IterPrint          = DebugPriLev > 0;
   snProb.PriLevOpt                   = DebugPriLev;
   snProb.GO.ProbL.optParam.IterPrint = DebugPriLev > 0;
   snProb.GO.ProbL.PriLevOpt          = DebugPriLev;
end

% HKH Use Prob.GO.ProbL for glcCluster and others, if LOCAL ~= 0

if(usecgolib == 1)
  gn_gFunc = 'gn_g';
else
  gn_gFunc = '';
end

gnProb = CGOGlobalProb(Prob, 'RBFgn', 'gn_f', gn_gFunc, [],...
   'rbf_c', 'rbf_dc', 'rbf_d2c', x_L, x_U, x_D, x_LL, x_UU, SCALE,...
   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
   backupSolver1, backupSolver2, PriSub);

gnProb.CGOLIB = Prob.CGOLIB;
gnProb.GO.ProbL.CGOLIB = Prob.CGOLIB;


% --------------------------------------------------
% Special input in gnProb
% --------------------------------------------------
if DebugPriLev > 0
   gnProb.optParam.IterPrint          = DebugPriLev > 0;
   gnProb.PriLevOpt                   = DebugPriLev;
   gnProb.GO.ProbL.optParam.IterPrint = DebugPriLev > 0;
   gnProb.GO.ProbL.PriLevOpt          = DebugPriLev;
end
gnProb.CGO.idea           = idea;
gnProb.GO.ProbL.CGO.idea  = idea;

% --------------------------------------------------

%--------------------------------------------------------------
% Define current best point
% (x_min, O_min, fMin) at index fIdx and flag Feasible
%--------------------------------------------------------------

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
if isempty(Cc)
   CcMin = [];
else
   CcMin = Cc(:,fIdx);
end

% Best point found in unit space
x_min = X(:,fIdx(1));

% Best point found in original space
if SCALE
   O_min = tomsol(9, x_L, x_min, x_D);
else
   O_min = x_min;
end
% --------------------------------------------------


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

OLDuc = 2;
if OLDuc >= 2
   snProbD1 = DefsnProbD1(snProb);
   snProbD1 = ProbCheck(snProbD1,localSolver,checkType('con'));
end
% Becase new sn_f,sn_g tests on CGO.fnStar
gnProb.CGO.fnStar = NaN;
if isempty(IntVars)
   snProb = ProbCheck(snProb,globalSolver,checkType('con'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('glc'));
else
   snProb = ProbCheck(snProb,globalSolver,checkType('minlp'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('minlp'));
end

if PriLev > 1
   fprintf('globalSolver %s ',globalSolver);
   fprintf('localSolver %s ',localSolver);
   fprintf('GOlocalSolver %s ',GOlocalSolver);
   fprintf('backupSolver1 %s ',backupSolver1);
   fprintf('backupSolver2 %s ',backupSolver2);
   fprintf('\n');
end

% ------------------------------------
% Lipschitz initialization
% ------------------------------------
% Find max and min Lipschitz constant in initial set

Dx     =  tomsol(30,X,X);
LipLow =  inf;
LipUpp = -inf;
for i=1:n-1
   L      = abs(F(i+1:n)-F(i))./Dx(i+1:n,i);
   LipLow = min(LipLow,min(L));
   LipUpp = max(LipUpp,max(L));
end


% *********************************************************
% *************** START MAIN ITERATION LOOP ***************
% *********************************************************

% n     = Number of points used in the interpolation
% nFunc = Total number of costly function evaluations
% nInit = Number of points used for startup, now or before
% Iter  = Number of fresh sampled points
% modN  = Cycle step

modN       = -1-max(0,infStep);
% addlocStep = 0; % NOT USED NOW
AddinfStep = 0;
LocRed     = 0;  % Keep track if local search points do reduce f(x)
GlobRed    = 0;  % Keep track if global search points do reduce f(x)
localGrid  = 0;  % Local Grid iteration number
globalGrid = 0;  % Global Grid iteration number

Iter       = 0;
Update     = 1;
fnStar     = Inf;
dSQRT      = sqrt(d); % Used in distance calculation

FLOWER     = fIdx(1);
fMinIter   = 0;
NOUPDATE   = 0;             % Used in test for Inform = 4
SUCCESS    = 1;             % Used to flag for No update possible in grid
SAME1      = 0;             % Used in test for Inform = 5
SAME2      = 0;             % Used in test for Inform = 6
MaxCycle   = 5;             % Used in tests for Inform 4,5,6
O_pre      = inf*ones(d,1);

% HKH New input parameters
% Gamma1 = distance between possible minima
Gamma1 = 0.03;

if isempty(xOptS)
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
   snOptg = norm(nlp_g(xOptS,snProb));
   if isnan(snOptg)
      snOptH = [];
      snOptE = NaN;
   else
      snOptH = nlp_H(xOptS,snProb);
      snOptE = sum(eig(snOptH)<0);
   end
   zz     = tomsol(30,xOptS,X);
   % Compute closest point to global optimum
   dXO    = min(zz);
   % Compute sum of all distances from global optimum to initial set X
   SumXO  = sum(zz);
   % Compute mean distance from global optimum to initial set X
   MeanXO = mean(zz);
end

Its = [];
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
   fprintf('\n');
   maxF = max(F);
   fprintf('  max(F) %9.5g', maxF);
   fprintf(' med(F) %9.5g', median(F));
   fprintf(' rng(F) %9.5g', maxF-fMin);
   fprintf(' pDist %9.5g ', pDist);
   fprintf(' LipUpp %f.',LipUpp)
   fprintf(' LipLow %f.',LipLow)
   fprintf('\n');
   xprint(O_min,'  xMin:',' %12.8f',8)
   if ~isempty(xOptS)
      if PriLev > 2
         xprint(Prob.x_opt,'  xOpt:',' %12.8f',8)
         % zz = tomsol(30,xOptS,X); % NOT USED NOW
         fprintf('  SumXO (sum||xOpt-X||) %f ',SumXO);
         fprintf('MeanXO (mean||xOpt-X||) %f ',MeanXO);
         doO = min(tomsol(30,x_min,xOptS));
         fprintf(' doO (||xOpt-xMin|| %f ',doO);
         fprintf('\n');
         if PriLev > 3
            fprintf('   dXO (min||xOpt-X||) %f ',dXO);
            fprintf('sn_f(xOpt) %8.5f ',snOptf);
            fprintf('||sn_g(xOpt)|| %f ',snOptg);
            fprintf('negeig(sn_H(xOpt)) %d ',snOptE);
            fprintf('\n');
         end
      end
   else
      Its.dXO    = [];
      Its.snOptf = [];
      Its.snOptg = [];
      Its.snOptE = [];
   end
end

% -------------- Result saving -------------------

TIME0     = Prob.TIME0;
DIGIT     = 0.1;
TESTDIG   = [];
convflag  = 0;
cpumax    = 0;
% progress  = 1;
% FLOW      = nFunc; % NOT USED NOW
LocErr    = Inf;
RelLocErr = Inf;
PROGRESS  = Inf*ones(N+1,1);

GlobErr(1:N)    = Inf;
RelGlobErr(1:N) = Inf;
% GlobImp(1:N)    = Inf; % NOT USED NOW

snProb.ixDist          = [];
snProb.GO.ProbL.ixDist = [];
gnProb.ixDist          = [];
gnProb.GO.ProbL.ixDist = [];
%HKH - needed?
vJFac = 1D0;

HKHPLOT = 0;
if HKHPLOT
   NamX = Name;
   NamX(NamX==' ') = '-';
   NamX(NamX==':') = '-';
   if DebugARBF
      fxFIG = figure('Name','fnStar and f(x(fnStar))');
      rbfT  = figure('Name','RBF Surface and x(fnS)');
      ProbT = plotProblem(Prob,O,F);
      set(ProbT,'Name','Function to Optimize and x(fnS)');
   end
   doMFIG1 = figure('Name','Distance1 minOfSurface to gnMin(Target value)');
   doMFIG2 = figure('Name','Distance2 minOfSurface to gnMin(Target value)');
   doMFIG3 = figure('Name','Distance3 minOfSurface to gnMin(Target value)');
end
%NHQ 
% PLOTFLAG 0/1. Flag used to control plots. 
% If set to 1, plots are produced when d == 2.
PLOTFLAG = 0;
if d ~= 2
   PLOTFLAG = 0;
end
if PLOTFLAG | (HKHPLOT & d==2)
   ProbFIG = plotProblem(Prob,O,F);
   set(ProbFIG,'Name','Function to Optimize');
   rbfFIG  = figure('Name','RBF Surface');
   gnFIG   = figure('Name','gnProb Function');
   gnProb.CGO.modN   = modN;
end

%NHQ 
% If RMSError == 1 compute RMS error of the interpolated surface.
RMSError = 0;

%NHQ  
% If ExDTEST == 1 compute ExD TEST Variables
ExDTEST  = 0;
minXdist = [];
xDIGITS  = [];
if ExDTEST == 1
   x_opt = Prob.x_opt';
   nrXopt = size(x_opt,2);
   %f_opt = Prob.f_opt(1);
   %x_box    = 0.1*x_D;
   x_box    = x_D/sqrt(10);

   x_opt_L = x_opt - repmat(0.5*x_box,1,nrXopt);
   x_opt_U = x_opt + repmat(0.5*x_box,1,nrXopt);

   % X_1% norm, any point inside box centered around x_opt?
   for k = 1:nrXopt
      while any( all(repmat(x_opt_L(:,k),1,nFunc) <= O) & ...
            all(repmat(x_opt_U(:,k),1,nFunc) >= O)  )
         xDIGITS = [xDIGITS [Iter nFunc]'];
         %x_box   = 0.1*x_box;
         x_box    = x_box/sqrt(10);
         x_opt_L = x_opt - repmat(0.5*x_box,1,nrXopt);
         x_opt_U = x_opt + repmat(0.5*x_box,1,nrXopt);
      end
   end
end

while nFunc < MaxFunc & convflag == 0

   % if nFunc-FLOW > MaxCycle*(N+1), progress = 0; break; end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end

   % Set parameters in global structure CGO
   gnProb.CGO.X          = X;
   gnProb.CGO.F          = Fpen;
   % gnProb.CGO.n        = n;
   % gnProb.CGO.d        = d;

   gnProb.CGO.rbfType    = rbfType;

   % Minimize s_n(y) using global optimization
   %  -  Solve  min[s_n(y)] with a local optimizer by starting from
   %     the interpolation point with the least function value

   Iter              = Iter + 1;


   % Determine fnStar value fnStarMax, limit where interpolation is WILD
   if fMin == 0
      v = min(F( F>1E-7 ));
      if isempty(v)
         v = 1E-7;
      end
      %HKH try
      % fnStarMax = -10*v;
      fnStarMax = -100*v;
   else
      %fnStarMax = fMin - 0.10*abs(fMin);
      fnStarMax = fMin - 0.50*abs(fMin);
   end
   % Set parameters used in snSolve
   snProb.snP              = Iter;
   snProb.x_0              = x_min;
   snProb.CGO              = gnProb.CGO;
   snProb.CGO.globalSolver = globalSolver;
   snProb.CGO.localSolver  = localSolver;
   %HKH NEW
   % modN not yet updated, make preliminary update for snSolve
   %snProb.CGO.modN         = mod(modN+1,N+1);
   snProb.CGO.modN         = modN;
   snProb.PriLev           = PriSub;

   %snProb.CGO.globalSolver = 'multiMin';
   %snProb.CGO.localSolver  = 'multiMin';
   snProb.CGO.fnStar       = NaN;
   snProb.GO.ProbL.CGO     = snProb.CGO;

   %NHQ plot!!!
   if PLOTFLAG
      % Update Sampled Points in figure ProbFIG
      figure(ProbFIG)
      plot3(O(1,end),O(2,end),F(end),'k.','MarkerSize',15)

      % Update and plot current RBF-surface
      %X = cgolib(109, DACEProb.CGOLIB.setid, n);
      snProb.probFile = Prob.probFile;
      clf(rbfFIG)
      plotProblem(snProb,X,[],rbfFIG);

      % Update and plot current gnProb
      %X = cgolib(109, DACEProb.CGOLIB.setid, n);
      modNtemp = gnProb.CGO.modN;
      gnProb.CGO.modN   = -1;
      gnProb.probFile = Prob.probFile;
      clf(gnFIG)
      plotProblem(gnProb,X,[],gnFIG);
      gnProb.CGO.modN   = modNtemp;
      pause(0.5)
      %keyboard
   end

   %if modN == -2 | modN == N
   %   RBFOK = 0;
   %else
   %   RBFOK = -2;
   %end
   RBFOK = 0;

   % FIX   =  0; % NOT MUCH USE OF FIX right now
   WILD  =  0;
   while RBFOK < 1
      snResult          = snSolve(snProb);
      min_sn            = snResult.f_k;
      if Feasible & min_sn-fMin > 1E-4*max(1,abs(fMin))
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
         if min_sn-fMin > 1E-4*max(1,abs(fMin))
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
      if PLOTFLAG
         disp('Min on sn-surface')
         disp(snResult.x_k)
         figure(rbfFIG)
         plot3(snResult.x_k(1),snResult.x_k(2),min_sn,'r*')
         pause(0.5)
         %keyboard
      end
      %NHQ min_sn < fnStarMax too easy??
      if min_sn < f_Low | min_sn < fnStarMax
         WILD = 1;
         if PriLev > 2
            fprintf('--- ENTER WILD MODE --- ');
            fprintf('fMin %12.6g ',fMin);
            fprintf('min_sn %12.6g ',min_sn);
            fprintf('fnStarMax %12.6g ',fnStarMax);
            fprintf('fLow %12.6g ',f_Low);
            fprintf('\n');
         end
      end

      if RBFOK == -2 & min_sn < f_Low
	if(usecgolib == 1)
	  error('Not implemented in cgolib');
	end
         %'TRY MEDIAN'
         tomsol(25) % Deallocates memory
         % F_S = F_m; % NOT USED
         F_m = min(median(F),F);
         control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
         % FIX   =  1;
         RBFOK = RBFOK+1;
      elseif RBFOK == -1 & min_sn < f_Low
	if(usecgolib == 1)
	  error('Not implemented in cgolib');
	end
         %'TRY LOG TRANSFORM'
         if fMin < 0
            F_m = log(1+fMin+F);
         else
            F_m = log(F);
         end
         tomsol(25) % Deallocates memory
         control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
         % FIX   =  1;
         RBFOK = RBFOK+1;
      else
         RBFOK = 1;
      end
   end
   if isempty(snResult.x_k)
      min_sn_y = x_min;
   else
      min_sn_y = snResult.x_k(:,1);
   end

   % c0 and min_snc NOT USED NOW
   %% Transform to original space
   %if SCALE
   %   min_snc = tomsol(9, x_L, min_sn_y, x_D);
   %   c0      = tomsol(9, x_L, x_min, x_D);
   %else
   %   min_snc = min_sn_y;
   %   c0      = x_min;
   %end


   % Choose a target value fnStar
   %if WILD & modN == -2 & fMinIter ~= Iter -1
   %   % Take infStep
   %   modN = -1;
   %elseif WILD
   %   modN = -2;
   %elseif modN == -2;
   %   modN = -infStep;  % If to start with infStep - maybe not good
   %   modN = 0;         % Always start with 1st step in cycle

   if modN == -8 | modN == -9 | modN == -10
      % Do nothing here, will use min_sn for this step
   elseif modN == -11 
      % Do nothing here, will take infStep before end of global grid with min_sn
   elseif modN == -5 | modN == -6 | modN == -7
      % Do nothing here, will use min_sn / Newton for this step
   elseif modN <= -3
      xLoc  = xSave;
   elseif WILD & modN == -2 & fMinIter ~= Iter-1
      % Take infStep
      modN = -2;
   elseif WILD
      modN = -2;
   elseif modN == -2;
      % modN = -infStep;  % If to start with infStep - maybe not good
      modN = 0;         % Always start with 1st step in cycle
   else
      modN = modN + 1;
   end
   if modN > N
      if FLOWER < nFunc - N
         % No improvement during the last cycle. Currently do nothing special
         modN = 0-max(0,infStep);
      else
         modN = 0-max(0,infStep);
      end
   end
   ucOK = 1;
   if modN == -1 | modN == -11
      fnStar = -inf;
   elseif modN == -2
      fnStar = fnStarMax;
   elseif modN == -8 | modN == -9 | modN == -10
      fnStar = fnStarMax;
   elseif modN == -5 | modN == -6 | modN == -7
      fnStar = fnStarMax;
   elseif modN == -3 | modN == -4
      fLoc   = fLocSave(1);
      fnStar = fnSSave(1);
   elseif modN == N
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
      else
         fnStar = min_sn;
      end
      gnProb.CGO.fnStar = fnStar;
   elseif  modN == 0 
      % GLOBAL GRID MODE - Global target value search
      SUCCESS = 1;
      globalGrid = globalGrid + 1;
      gnProb.CGO.modN   = modN;
      gnProb.x_0        = x_min;
      % Set iteration dependent parameters used in gnSolve, output from snSolve
      gnProb.IX   = snResult.multiMin.IX;
      gnProb.fOpt = snResult.multiMin.fOpt;
      gnProb.xOpt = snResult.multiMin.xOpt;
      vJ  = [0,1E-7,1E-6,1E-5,5E-5,1E-4,5E-4,1E-3,3E-3,5E-3,7E-3,8.5E-3,0.01*[1:9],0.1*[1:9],[1:9],...
             10,15,20,30,40,50,65,80,100,200,500,1000,5000,10000, inf];
      vJFac = max(1E-8,abs(min_sn));
      fSpan   = max(F_m)-fMin;
      if fMin > 0
         fRange = min(max(fMin,1),fSpan);
      else
         fRange = min(10*max(abs(fMin),1),fSpan);
      end
      if PriLev > 2
         fprintf('\n');
         fprintf('--- GLOBAL ---  TARGET VALUE SEARCH');
         fprintf('\n    ');
         fprintf('fMin       %12.6g ',fMin);
         fprintf('min_sn     %12.6g ',min_sn);
         fprintf('fnStarMax  %12.6g ',fnStarMax);
         fprintf('fSpan      %12.6g ',fSpan);
         fprintf('fRange     %12.6g ',fRange);
         fprintf('vJFac      %8.3g ',vJFac);
         fprintf('\n');
         xprinti(min_sn_y,'    min_sn_y:');
      end
      % Setup for search in target values
      maxfnS  = 90;
      XJ      = zeros(d,maxfnS);
      if isempty(IntVars)
         nfJ     = 8+DebugARBF;
      else
         nfJ     = 11+DebugARBF;
      end
      fJ      = zeros(maxfnS,nfJ);
      fJ(:,4) = NaN; % Mark every row as NOT COMPUTED

      % Columns in fJ (fV row), quantities computed for xTry
      % 1: fnStar 2: onB(xTry) 3: s_n(xTry) 4: gn_f(xTry) 5: my(xTry)
      % 6: doX 7: doM 8: doS if ~isempty(IntVars): 9: doX(2) 10: doM(2) 11: doS(2)
      % Grid of target values
      fnS                = min_sn-vJ*vJFac;
      mTry               = length(fnS);
      % First value is minimum of RBF surface
      xTry               = min_sn_y;
      fTry               = min_sn;
      XJ(:,1)            = xTry;
      fV                 = fnSStat(fnS(1),xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
      fJ(1,1:nfJ-DP)     = fV;
      if DebugARBF
         fx              = fTrue(xTry,x_L,x_D,Prob,varargin{:});
         fJ(1,nfJ)       = fx;
      end
      minonB             = fV(2);
      % Compute last value - fnStar = -inf
      gnProb.CGO.fnStar  = 0;
      gnProb.CGO.alpha   = 0;
      gnProb.CGO.modN    = -1;
      gnProb.x_0         = min_sn_y;
      gnProb.X0          = XJ(:,1);
      gnProb.GO.ProbL.CGO = gnProb.CGO;
      gnRInf             = gnSolve(gnProb);
      gnProb.CGO.modN    = modN;
      xTry               = gnRInf.x_k(:,1);
      fTry               = gnRInf.f_k;
      XJ(:,maxfnS)       = xTry;
      fV                 = fnSStat(fnS(end),xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
      fJ(maxfnS,1:nfJ-DP)= fV;
      if DebugARBF
         fx              = fTrue(xTry,x_L,x_D,Prob,varargin{:});
         fJ(maxfnS,nfJ)  = fx;
      end
      doMInf             = fJ(maxfnS,7);
      fTryInf            = fJ(maxfnS,4);
      XX                 = XJ;
      if ~SCALE
         XX(Reals,1)      = (XJ(Reals,1)-x_L(Reals))./x_D(Reals);
         XX(Reals,maxfnS) = (XJ(Reals,maxfnS)-x_L(Reals))./x_D(Reals);
      end
      FASTSTOP           = inf;
      SetvJ              = 0;
      fnStarJump         = inf;
      nTry               = 1;
      j                  = 1;
      while j < mTry-1
         j                 = j+1;
         fnStar            = fnS(j);
         if isnan(fJ(j,4))
            % No solution is computed for the current fnS(j)
            gnProb.CGO.fnStar = fnStar;
            % Use minimum on surface as initial x_0 for current target value
            gnProb.x_0        = min_sn_y;
            % Use previous solutions for all target values as matrix of initial points
            gnProb.X0         = XJ(:,[maxfnS,1:nTry]);
            gnProb.GO.ProbL.CGO = gnProb.CGO;
            gnR0              = gnSolve(gnProb);
            xTry              = gnR0.x_k(:,1);
            fTry              = gnR0.f_k;
            fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
            % Compute solution scaled to [0,1]
            XXTry            = xTry;
            if ~SCALE
               XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
            end
            % Accept new point and update
            nTry              = nTry+1;
            XJ(:,nTry)        = xTry;
            XX(:,nTry)        = XXTry;
            fJ(nTry,1:nfJ-DP) = fV;
            if DebugARBF
               fx             = fTrue(xTry,x_L,x_D,Prob,varargin{:});
               fJ(nTry,nfJ)   = fx;
            end
         else
            % A solution was previously computed for the current fnS(j)
            nTry              = nTry+1;
            fV                = fJ(nTry,1:nfJ-DP);
            xTry              = XJ(:,nTry);
            fTry              = fV(4);              % Same as fJ(nTry,4)
            XXTry             = XX(:,nTry);
         end
         doMNew            = fV(7);
         % Check if previous solutions was OK
         % needed? HKH
         for k=nTry-1:2
            doMOld = fJ(k,7);
            if doMOld - doMNew > 1E-15 * doMNew
               if PriLev > 1
                  fprintf('    Must recompute solution %d\n',k);
               end
               % Use solution from new target value as initial point
               gnProb.x_0        = xTry;
               % Use best minimum in previous computation as additional initial value
               gnProb.X0         = XJ(:,k);
               % HKH0804 - must set correct fnStar
               gnProb.CGO.fnStar = fnS(k);
               gnProb.GO.ProbL.CGO = gnProb.CGO;

               gnR0              = gnSolve(gnProb);
               xTry              = gnR0.x_k(:,1);
               fTry              = gnR0.f_k;
               fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
               % Compute solution scaled to [0,1]
               XXTry            = xTry;
               if ~SCALE
                  XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
               end
               XJ(:,k)        = xTry;
               XX(:,k)        = XXTry;
               fJ(k,1:nfJ-DP) = fV;
               % Should also print out something
            end
         end
         D                 = normDist(XXTry,XX(:,nTry-1),dSQRT,IntVars,Reals);
         onB               = fV(2);
         if D(1) > 3*Gamma1 & onB <= minonB & onB == fJ(nTry-1,2)
            % Big difference between solutions, try to improve solution in case of bad optimum
            if PriLev > 1
               xprinti(fJ(1:nTry,2),'onB:');
               fprintf('    Recompute solution to possibly decrease D %f',D(1));
            end
            % Use solution from new target value as initial point
            gnProb.x_0        = xTry;
            % Use best minimum in previous computation as additional initial value
            gnProb.X0         = XJ(:,nTry-1);
            gnProb.CGO.fnStar = fnStar;
            gnProb.GO.ProbL.CGO = gnProb.CGO;
            gnR0              = gnSolve(gnProb);
            xTry              = gnR0.x_k(:,1);
            fTry              = gnR0.f_k;
            fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
            % Compute solution scaled to [0,1]
            XXTry             = xTry;
            if ~SCALE
               XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
            end
            XJ(:,nTry)        = xTry;
            XX(:,nTry)        = XXTry;
            fJ(nTry,1:nfJ-DP) = fV;
            D = normDist(XXTry,XX(:,nTry-1),dSQRT,IntVars,Reals);
            onB               = fV(2);
            if PriLev > 1
               % Should also print out something
               fprintf(' --- D now %f\n',D(1));
            end
         end
         doMOld            = fJ(nTry-1,7);
         fTryOld           = fJ(nTry-1,4);
         if PriLev > 1 & doMNew - doMInf > 1E-5 * doMInf
            fprintf('%2d: fnStar %15.7f doMInf %20.15f < doMNew %20.15f\n',j,fnStar,doMInf, doMNew);
            fprintf('  :               fTryInf %25.18f   fTryNew %25.18f\n',fTryInf,fTry);
    %HKHNY Possibly recompute inf solution, check if f-vales are same 1st
         end
         if SetvJ & (doMNew - doMInf >= 1E-15 * doMInf | onB == d)
            if PriLev > 1
               fprintf('    ');
               fprintf('SetvJ: %2d: fNSTAR %15.7f doMInf %20.15f < doMNew %20.15f',j,fnStar,doMInf, doMNew);
               fprintf(' Decrease vJFac from %12.6f',vJFac);
            end
            vJVal = (min_sn - fnStar) / (vJFac*fRange);
            if RelLocErr < 1 & vJFac > 1E-3
               % if vJVal <= 1, vJFac =vJVal/vJ(end-2) * vJFac; end
               if vJVal <= 1, vJFac =max(vJFac/10,vJVal/vJ(end-2) * vJFac); end
               if PriLev > 1
                  fprintf(' to %12.6f',vJFac);
               end
            else
               if PriLev > 1
                  fprintf(' ! No Decrease, RelLocErr %f',RelLocErr);
               end
            end
            SetvJ = 0;
            if PriLev > 1
               fprintf('\n');
            end
         end
         if SetvJ & (nTry >= mTry-1)
            if vJFac < 1E3 & onB <= minonB
               vJFac = vJFac*10;
               SetvJ = 0;
               if PriLev > 1
                  fprintf('%2d: Increase vJFac 10 times to %20.15f\n',j,vJFac);
               end
             end
         end
         if PriLev > 1 & doMNew - doMInf > 1E-15 * doMInf
            fprintf('%2d: fNSTAR %15.7f doMInf %20.15f < doMNew %20.15f\n',j,fnStar,doMInf, doMNew);
    %HKHNY Possibly recompute inf solution, check if f-vales are same 1st
    %HKHNY 2 sådana tester, varför
         end
         %if doMOld - doMNew > 1E-15 * doMNew
         %   fprintf('%2d: fNSTAR %15.7f doMNew %20.15f < doMOld %20.15f\n',j,fnStar,doMNew, doMOld);
         %end
         if PriLev > 1 & doMOld - doMNew > 1E-5 * doMNew
            fprintf('%2d: fnStar %15.7f doMNew %20.15f < doMOld %20.15f\n',j,fnStar,doMNew, doMOld);
            fprintf('  :               fTryNew %25.18f   fTryOld %25.18f\n',fTry,fTryOld);
    %HKHNY Kolla om fvalues lika eller ej
         end
    %HKHNY Skriv ut om flera min med samma f-vals, o deras doM
         M = 0;
         if D(1) > 0 & mTry < maxfnS-1 & fnStar < fnStarJump
            if PriLev > 1
               fprintf('%2d: D %12.6f [%d]',j,D(1), onB)
               if D(1) > Gamma1 & onB == fJ(nTry-1,2), fprintf('    Too big range'); end
               fprintf('\n')
            end
            if D(1) > Gamma1 & onB <= minonB & onB == fJ(nTry-1,2)
               %HKH M = ceil(D(1)/Gamma1);
               M = 1+ceil(D(1)/Gamma1);
               if PriLev > 1
                  fprintf('    Must divide further, %d more points\n',M);
               end
               %HKH M = min([5,maxfnS-mTry-1,M]);
               M = min([7,maxfnS-mTry-1,M]);
               fnStarJump = fnStar;
            end
         end
         if onB > fJ(nTry-1,2) & D(1) > Gamma1 & fnStar < fnStarJump
            % New solution has more components on bounds
            if onB == fJ(nTry-1,2)+1 & D(1) > Gamma1 & onB < fJ(maxfnS,2)
               %HKH M = ceil(D(1)/Gamma1);
               M = 1+ceil(D(1)/Gamma1);
               if PriLev > 1
                  fprintf('    Must divide further II, %d more points\n',M);
               end
               %HKH M = min([5,maxfnS-mTry-1,M]);
               M = min([7,maxfnS-mTry-1,M]);
               fnStarJump = fnStar;
            elseif onB == fJ(nTry-1,2)+1
               % Normal change for new solution
            else
               % More than one component jumped to a bound
               if PriLev > 1
                  fprintf('    More than one component jumped to a bound: %d to %d \n',fJ(nTry-1,2),onB);
               end
               if nTry <= 10
                  M = min([10,maxfnS-mTry-1,3*(onB-fJ(nTry-1,2))]);
               else
                  M = min([5,maxfnS-mTry-1,onB-fJ(nTry-1,2)]);
               end
               fnStarJump = fnStar;
            end
         end
         if M > 0
            fnStarRange = fnS(nTry-1)-fnStar;
            Newfn = fnS(nTry-1)-[1:M]/(M+1)*fnStarRange;
            if PriLev > 1
               fprintf('    Set M = %d more points',M);
               fprintf(' fnStarRange %f',fnStarRange);
               fprintf(' fnStar %f',fnStar);
               fprintf(' fnSnTry %f',fnS(nTry-1));
               fprintf('\n');
               %Newfn = fnS(nTry-1)-log([1:M])/log((M+1))*fnStarRange;
               xprint(Newfn,'    Newfn:');
            end
            % Add new fnStar values first in list, move current to new position
            fnS  = [fnS(1:j-1),Newfn,fnS(j:end)];
            vJ   = [vJ(1:j-1),ones(1,M)*vJ(j),vJ(j:end)];
            % Move current fnStar solution to new position
            XJ(:,nTry+M)        = XJ(:,nTry);
            XX(:,nTry+M)        = XX(:,nTry);
            fJ(nTry+M,:)        = fJ(nTry,:);
            % Mark row nTry as not computed
            fJ(nTry,4)          = NaN;
            mTry                = mTry + M;
            nTry                = nTry - 1;
            j                   = j - 1;
         end
         % Avoid that range of fnStar values down to -inf with interesting solutions is not utilized
         if nTry == mTry-1 & onB ~= fJ(maxfnS,2)
            M                = maxfnS-mTry;
            mTry             = maxfnS;
            fnStarRange1     = fnS(1)-fnStar;
            fnStarRange2     = fnS(nTry-1)-fnStar;
            fnStarRange      = max(3*2^(-M)*fnStarRange1,fnStarRange2);
            Newfn            = fnS(nTry)-2.^[0:M-1]*fnStarRange;
            fprintf('    Increase Range to -inf:');
            fprintf(' fnStar %f',fnStar);
            fprintf(' fnStarRange1 %f',fnStarRange1);
            fprintf(' fnStarRange2 %f',fnStarRange2);
            fprintf(' fnStarRange %f',fnStarRange);
            fprintf('\n');
            xprint(Newfn,'    Newfn:');
            % Add new fnStar values between current and -inf in list
            fnS  = [fnS(1:j),Newfn,fnS(end)];
            vJ   = [vJ(1:j),ones(1,M)*vJ(j),vJ(end)];
         end
         %if nTry == mTry-1 & onB == fJ(maxfnS,2) & onB == minonB
         if nTry == mTry-1 & onB == fJ(maxfnS,2)
            D  = normDist(XXTry,XX(:,maxfnS),dSQRT,IntVars,Reals);
            if D(1) > Gamma1
               fprintf('    Here we should add some items down to -inf. D= %f',D(1));
               fprintf('\n');
            end
         end
         if FASTSTOP == nTry, break; end
         if FASTCODE & isinf(FASTSTOP) & ...
            (onB <= fJ(maxfnS,2) & abs(doMNew - doMInf) <= 1E-5 * doMInf )
            %(onB <= fJ(maxfnS,2) & onB > minonB & abs(doMNew - doMInf) <= 1E-5 * doMInf )
            %((onB == fJ(maxfnS,2) & onB > minonB) | (doMNew - doMInf > 1E-5 * doMInf & onB < fJ(maxfnS,2)))
            if M==0
               break;
            else
               FASTSTOP = nTry+M+1;
            end
         end
      end
      % Move -inf row to row nTry+1
      if nTry ~= maxfnS-1
         nTry              = nTry+1;
         XJ(:,nTry)        = XJ(:,maxfnS);
         XX(:,nTry)        = XX(:,maxfnS);
         fJ(nTry,:)        = fJ(maxfnS,:);
      end
      % Compute x_min solution scaled to [0,1]
      x_minS = x_min;
      if ~SCALE
         x_minS(Reals) = (x_minS(Reals)-x_L(Reals))./x_D(Reals);
      end
      D0                = normDist(x_minS,XX(:,1),dSQRT,IntVars,Reals);
      if PriLev > 1
         fprintf('Scaled distance from x_min to 1st solution (min_sn_y): D0 %f\n',D0(1));
      end

      % Print one row for each fnStar tried
      if PriLev > 2
         fnSPrint(XJ(:,1:nTry),fJ(1:nTry,:),min_sn,DebugARBF)
      end

      onBTry    = fJ(1:nTry,2);
      [g,DistGrp] = djCluster(XX(:,1:nTry),onBTry,epsX,PriLev > 2,IntVars,Reals);
      nGroup = g(end);

      fTry      = fJ(1:nTry,4);
      if DebugARBF
         fxTry      = fJ(1:nTry,end);
      else
         fxTry      = [];
      end
      doXTry    = fJ(1:nTry,6);
      % doMTry    = fJ(1:nTry,7);

      % xTry    = zeros(nGroup,1); % SEEMS NOT USED
      Gsize   = zeros(nGroup,1);
      GonBmin = zeros(nGroup,1);
      GonBmax = zeros(nGroup,1);
      ii      = 0;
      while ii < nGroup
         ii = ii+1;
         ix = find(g==ii);
         Gsize(ii) = length(ix);
         minonB      = min(onBTry(ix));
         GonBmin(ii) = minonB;
         GonBmax(ii) = max(onBTry(ix));
      end
      if PriLev > 1
         fprintf('RelGlobErr-1 %f ',RelGlobErr(1));
         fprintf('GlobErr-1 %f ',GlobErr(1));
         fprintf('\n');
         fprintf('RelGlobErr-2 %f ',RelGlobErr(2));
         fprintf('GlobErr-2 %f ',GlobErr(2));
         fprintf('\n');
         fprintf('RelLocErr    %f ',RelLocErr);
         fprintf('LocErr    %f ',LocErr);
         fprintf('\n');
         fprintf('\n');
      end
      xCand   = [];
      minonBG = min(GonBmin);
      % maxonBG = min(GonBmax); % NOT USED
      for i=1:nGroup
         j  = i;
         ix = find(g==j);
         iy = find(onBTry(ix) == GonBmin(j));
         iz = ix(iy);
         minDist  = min(DistGrp(iz));
         L  = max(DistGrp(iz))- minDist;
         % L  = max(DistGrp(ix));
         % l  = ceil(10*L); % HKH NOTE l SEEMS NOT TO BE USED
         if j == 1
            M =  max(2,ceil(L/Gamma1));
            if D0(1) > L | D0(1) >= min(RelLocErr,Gamma1/3)
               D01=D0(1);
               k       = iz(1);
               xCand = [xCand,k];
               if PriLev > 1
                  fprintf('-------- WHY1? D01 %f ',D01)
                  fprintf('L %d ',L)
                  fprintf('M %d ',M)
                  fprintf('RelLocErr %e',RelLocErr)
                  fprintf('\n')
                  GrpPrint(i,j,k,'xCand1-DO',GonBmin,GonBmax,Gsize,fxTry,doXTry)
               end
            end
            % Generate M points, the points {ii/M*GrpDistance, i=1,M}
            for ii = 1:M-1
               %[vv,iy] = min(abs(DistGrp(iz)-(minDist+ii*L/M)));
               [vv,iy] = min(abs(DistGrp(iz)-(ii*L/M)));
               k       = iz(iy);
               if any(xCand==k) & length(iz) > iy
                  k       = iz(iy+1);
               end
               if ~any(xCand==k)
                  xCand = [xCand,k];
                  if PriLev > 1
                     GrpPrint(i,j,k,['xCand1-D',num2str(ii)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
                  end
               end
               %k = iz(ii);
               %if DistGrp(k)-D > 0.03
               %   xCand   = [xCand,k];
               %   D       = DistGrp(k);
               %   GrpPrint(i,j,k,'xCandmZ',GonBmin,GonBmax,Gsize,fxTry,doXTry)
               %end
            end
            k       = iz(end);
            xCand = [xCand,k];
            if PriLev > 1
               GrpPrint(i,j,k,['xCand1-D',num2str(M)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
            end
            %elseif j == nGroup
         elseif (j < nGroup & minonBG == GonBmin(i)) | (j == nGroup & minonBG == GonBmin(i))
            % Right now both inbetween groups 2,...,nGroup-1, as well as last group
            % elseif minonBG == GonBmin(i)
            M =  max(1,ceil(L/Gamma1));
            if M > 2
               k       = iz(1);
               xCand = [xCand,k];
               if PriLev > 1
                  GrpPrint(i,j,k,['xCand',num2str(i),'-1'],GonBmin,GonBmax,Gsize,fxTry,doXTry)
               end
            end
            % Generate M points, the points {ii/M*GrpDistance, i=1,M}
            for ii = 1:M-1
               %[vv,iy] = min(abs(DistGrp(iz)-(minDist+ii*L/M)));
               [vv,iy] = min(abs(DistGrp(iz)-(ii*L/M)));
               k       = iz(iy);
               if any(xCand==k) & length(iz) > iy
                  k       = iz(iy+1);
               end
               if ~any(xCand==k)
                  xCand = [xCand,k];
                  if PriLev > 1
                     GrpPrint(i,j,k,['xCand',num2str(i),'-',num2str(ii)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
                  end
               end
            end
            k       = iz(end);
            xCand = [xCand,k];
            if PriLev > 1
               GrpPrint(i,j,k,['xCand',num2str(i),'-',num2str(M)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
            end
         elseif 1 & j == nGroup
            % HKH US Special
            % Take 1st item in last group
            k     = iz(1);
            xCand = [xCand,k];
            if PriLev > 1
               GrpPrint(i,j,k,['xCand',num2str(i),'-inf'],GonBmin,GonBmax,Gsize,fxTry,doXTry)
            end
         end
      end
      if 0 % HKHUS SPECIAL ADD
         xCand = [xCand,nTry];
      end
      ix = find(doXTry(xCand) > epsX);
      if isempty(ix)
         xCand = xCand(length(xCand));
         nCand = 1;
      elseif length(ix) ~= length(xCand)
         nCand  = length(xCand);
         if PriLev > 1
            fprintf('Remove %d candidates too close to X\n',nCand-length(ix))
         end
         xCand = xCand(ix);
         nCand  = length(xCand);
      else
         nCand  = length(xCand);
      end

      xLoc   = XJ(:,xCand);
      if nCand > 1
         if PriLev > 1
            xprinti(xCand,'xCand');
         end
         xC = xCand;
         for i = nCand:-1:2
            ii = xCand(i);
            for j=1:i-1
               jj = xCand(j);
               if tomsol(30,XX(:,jj),XX(:,ii))/dSQRT < 1E-6
                  xLoc = [xLoc(:,1:i-1),xLoc(:,i+1:end)];
                  xC   = [xC(1:i-1),xC(i+1:end)];
                  if PriLev > 1
                     fprintf('    Remove candidate %d, #%d, is too close to #%d\n',i,ii,jj);
                  end
                  break;
               end
            end
         end
         if PriLev > 1 & length(xC) < nCand
            xprinti(xC,'    Compressed xCand');
         end
         xCand    = xC;
         nCand    = size(xLoc,2);
         fLocSave = fTry(xCand(1:end));
         fnSSave  = fnS(xCand(1:end));
      end

      fLoc    = fTry(xCand(1));
      fnStar  = fnS(xCand(1));
      GlobRed = 0;    % Keep track if global points do reduce f(x)
      LocRed  = NaN;  % No local grid now
      if nCand > 1
         if PriLev > 1
            fprintf('    Try %d new points. ',nCand);
         end
         modN = -3;
      elseif nCand == 1
         if PriLev > 1
            fprintf('    Try 1 new point.  ');
         end
         modN = -3;
      else
         if PriLev > 1
            fprintf('    No global points found.  ');
         end
         fnStar = min_sn;
         modN   = -8;
      end
      if PriLev > 1
         fprintf('\n');
         if nCand > 0
            for i = 1:nCand
               %fprintf('xNew #%d:',xCand(i));
               %fprintf('\n');
               xprint(xLoc(:,i),'  xNew:');
            end
         end
      end
      %HKH plot!!!
      if HKHPLOT
         ixG = find(fJ(1:nTry,2) == 0);
         ixR = find(fJ(1:nTry,2) == d);
         ixB = find(fJ(1:nTry,2) < d & fJ(1:nTry,2) > 0);
         if 1 & DebugARBF
            x1 = 1:nTry; %fnStar
            y1 = fJ(1:nTry,end); %fx
            y2 = fJ(1:nTry,1);   %fnStar
            clf(fxFIG)
            figure(fxFIG)
            plot(x1,y2,'c-')
            plot(x1,y2,'cx')
            hold on
            plot(x1,y1,'k-')
            if ~isempty(ixG)
               plot(x1(ixG),y1(ixG),'go')
            end
            if ~isempty(ixR)
               plot(x1(ixR),y1(ixR),'ro')
            end
            if ~isempty(ixB)
               plot(x1(ixB),y1(ixB),'bo')
            end
            ylabel('f_n^* (cyan) and f(x) (black - green blue red)')
            title(['f(x) and f_n^*.  Global Iter ' num2str(globalGrid),...
                   ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
            print( fxFIG, '-djpeg', [NamX,'fxGoIt',num2str(globalGrid),...
	                               'x' num2str(rbfType)]); 
            hold off
         end
         x1 = fJ(1:nTry,1); %fnStar
         y1 = fJ(1:nTry,7); %doM
         y2 = fJ(1:nTry,6); %doX

if 1
         % Plot doM
         clf(doMFIG1)
         figure(doMFIG1)
         % Update Sampled Points in figure ProbFIG
         plot(x1,y1,'k-')
         %plot(x1,y1,'cx')
         hold on
         %plot(x1,y2,'b-')
         if ~isempty(ixG)
            plot(x1(ixG),y1(ixG),'go')
         end
         if ~isempty(ixR)
            plot(x1(ixR),y1(ixR),'ro')
         end
         if ~isempty(ixB)
            plot(x1(ixB),y1(ixB),'bo')
         end
         xlabel('f_n^*')
         ylabel('distance || x_{Min} - x(f_n^*) ||')
         title(['Distance x_{Min} to x(f_n^*). Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         print( doMFIG1, '-djpeg', [NamX,'doM1GoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 
end

if 1
         x1=fMin+1E-10-fJ(1:nTry,1);
if any(x1 <0),
disp(x1)
x1 = max(0,x1);
%keyboard
end
         clf(doMFIG2)
         figure(doMFIG2)
         % Update Sampled Points in figure ProbFIG
         semilogx(x1,y1,'k-')
         %semilogx(x1,y1,'cx')
         hold on
         %semilogx(x1,y2,'b-')
         if ~isempty(ixG)
            semilogx(x1(ixG),y1(ixG),'go')
         end
         if ~isempty(ixR)
            semilogx(x1(ixR),y1(ixR),'ro')
         end
         if ~isempty(ixB)
            semilogx(x1(ixB),y1(ixB),'bo')
         end
         xlabel('log10(f_{Min} - f_n^*)')
         ylabel('distance || x_{Min} - x(f_n^*) ||')
         title(['Distance x_{Min} to x(f_n^*). Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         print( doMFIG2, '-djpeg', [NamX,'doM2GoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 
         hold off
end

if 1
if any(y1 <0),
disp(y1)
y1=max(y1,0);
%keyboard
end
         clf(doMFIG3)
         figure(doMFIG3)
         loglog(x1,y1,'k-')
         %loglog(x1,y1,'cx')
         hold on
         %loglog(x1,y2,'b-')
         if ~isempty(ixG)
            loglog(x1(ixG),y1(ixG),'go')
         end
         if ~isempty(ixR)
            loglog(x1(ixR),y1(ixR),'ro')
         end
         if ~isempty(ixB)
            loglog(x1(ixB),y1(ixB),'bo')
         end
         xlabel('log10(f_{Min} - f_n^*)')
         ylabel('log10distance || x_{Min} - x(f_n^*) ||')
         title(['Distance x_{Min} to x(f_n^*). Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         print( doMFIG3, '-djpeg', [NamX,'doM3GoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 
         hold off
end

%if 1
if d==2
         snProb.probFile = Prob.probFile;
         clf(rbfT)
         plotProblem(snProb,XX,[],rbfT);
         hold on
         xlabel('x_1(f_n^*)')
         ylabel('x_2(f_n^*)')
         title(['x(f_n^*) on RBF. Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         %hold off
         print( rbfT, '-djpeg', [NamX,'rbfTGoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 

         figure(ProbT)
         %plot3(O(1,end),O(2,end),F(end),'k.','MarkerSize',15)
         z1 = fJ(1:nTry,end); %fx
         plot3(XJ(1,1:nTry),XJ(2,1:nTry),z1,'k.','MarkerSize',15)
         hold on
         xlabel('x_1(f_n^*)')
         ylabel('x_2(f_n^*)')
         title(['x(f_n^*) on f(x) surface. Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         %hold off
         print( ProbT, '-djpeg', [NamX,'probTGoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 

         figure(ProbFIG)
         plot3(O(1,end),O(2,end),F(end),'k.','MarkerSize',15)
         hold on
         xlabel('x_1')
         ylabel('x_2')
         title(['f(x) surface. Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         %hold off
         print( ProbFIG, '-djpeg', [NamX,'probGoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 

         % Update and plot current RBF-surface
         %X = cgolib(109, DACEProb.CGOLIB.setid, n);
         snProb.probFile = Prob.probFile;
         clf(rbfFIG)
         plotProblem(snProb,X,[],rbfFIG);
         hold on
         xlabel('x_1')
         ylabel('x_2')
         title(['RBF surface. Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         %hold off
         print( rbfFIG, '-djpeg', [NamX,'rbfGoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 

         % Update and plot current gnProb
         %X = cgolib(109, DACEProb.CGOLIB.setid, n);
         modNtemp = gnProb.CGO.modN;
         gnProb.CGO.modN   = -1;
         gnProb.probFile = Prob.probFile;
         clf(gnFIG)
         plotProblem(gnProb,X,[],gnFIG);
         gnProb.CGO.modN   = modNtemp;
         hold on
         xlabel('x_1')
         ylabel('x_2')
         title(['f_n^*= -\infty surface. Global Iter ' num2str(globalGrid),...
                ' Iter ' num2str(Iter) ' FuncEv ',num2str(nFunc)]);
         hold off
         print( gnFIG, '-djpeg', [NamX,'gnGoIt',num2str(globalGrid),...
	                            'x' num2str(rbfType)]); 
         pause(0.5)
         %keyboard
end
      end
%if globalGrid == 2 
%keyboard
%end
      AddinfStep = infStep < 0;

   elseif modN == 1

      localGrid = localGrid + 1;
      SUCCESS = 1;
      % LOCAL GRID MODE
      % Local target value search
      gnProb.CGO.modN   = modN;
      gnProb.x_0        = x_min;
      % Set iteration dependent parameters used in gnSolve, output from snSolve
      gnProb.IX   = snResult.multiMin.IX;
      gnProb.fOpt = snResult.multiMin.fOpt;
      gnProb.xOpt = snResult.multiMin.xOpt;
      % vJ  = [0 1E-4 1E-3, 0.01*[1:13],0.15,0.20,0.25 0.30 0.40 0.50 0.75 1.00 1.5 2 3 1E2 -inf];
      %vJ  = [0 1E-4 2.5E-4 5E-4 7.5E-4 1E-3, 2.5E-3 5E-3 7.5E-3, 0.01*[1:13],0.15,0.20,0.25 0.30 0.40 0.50 0.75 1.00 1.5 2 3 1E2 -inf];
      %vJ  = [0,1E-7,1E-6,1E-5,5E-5,1E-4,5E-4,1E-3,3E-3,7E-3,0.01*[1:9],0.1*[1:9],[1:9],...
      %       10,20,40,70,100,500,1000,5000,10000, inf];
      vJ  = [0,1E-7,5E-7,1E-6,5E-6,1E-5,3E-5,5E-5,7E-5,1E-4*[1:9],1E-3*[1:9],0.01*[1:9],0.1*[1:9],[1:2:9],...
             10,50,100,inf];
      vJFac = max(1E-8,abs(min_sn));
      fSpan   = max(F_m)-fMin;
      if fMin > 0
         fRange = min(max(fMin,1),fSpan);
      else
         fRange = min(10*max(abs(fMin),1),fSpan);
      end
      if PriLev > 2
         fprintf('\n');
         fprintf('---- LOCAL ---  TARGET VALUE SEARCH');
         fprintf('\n    ');
         fprintf('fMin       %12.6g ',fMin);
         fprintf('min_sn     %12.6g ',min_sn);
         fprintf('fnStarMax  %12.6g ',fnStarMax);
         fprintf('fSpan      %12.6g ',fSpan);
         fprintf('fRange     %12.6g ',fRange);
         fprintf('vJFac      %8.3g ',vJFac);
         fprintf('\n');
         xprinti(min_sn_y,'    min_sn_y:');
      end
      % Setup for search in target values
      maxfnS  = 90;
      XJ      = zeros(d,maxfnS);
      if isempty(IntVars)
         nfJ     = 8+DebugARBF;
      else
         nfJ     = 11+DebugARBF;
      end
      fJ      = zeros(maxfnS,nfJ);
      fJ(:,4) = NaN; % Mark every row as NOT COMPUTED

      % Columns in fJ (fV row), quantities computed for xTry
      % 1: fnStar 2: onB(xTry) 3: s_n(xTry) 4: gn_f(xTry) 5: my(xTry)
      % 6: doX 7: doM 8: doS if ~isempty(IntVars): 9: doX(2) 10: doM(2) 11: doS(2)

      % Grid of target values
      fnS                = min_sn-vJ*vJFac;
      mTry               = length(fnS);
      % First value is minimum of RBF surface
      xTry               = min_sn_y;
      fTry               = min_sn;
      XJ(:,1)            = xTry;

      fV                 = fnSStat(fnS(1),xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
      fJ(1,1:nfJ-DP)     = fV;
      if DebugARBF
         fx              = fTrue(xTry,x_L,x_D,Prob,varargin{:});
         fJ(1,nfJ)       = fx;
      end
      minonB             = fV(2);
      % Compute last value - fnStar = -inf
      gnProb.CGO.fnStar  = 0;
      gnProb.CGO.alpha   = 0;
      gnProb.CGO.modN    = -1;
      gnProb.x_0         = min_sn_y;
      gnProb.X0          = XJ(:,1);
      gnProb.GO.ProbL.CGO = gnProb.CGO;
      gnRInf             = gnSolve(gnProb);
      gnProb.CGO.modN    = modN;
      xTry               = gnRInf.x_k(:,1);
      fTry               = gnRInf.f_k;
      XJ(:,maxfnS)       = xTry;
      fV                 = fnSStat(fnS(end),xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
      fJ(maxfnS,1:nfJ-DP)= fV;
      if DebugARBF
         fx              = fTrue(xTry,x_L,x_D,Prob,varargin{:});
         fJ(maxfnS,nfJ)  = fx;
      end
      doMInf             = fJ(maxfnS,7);
      fTryInf            = fJ(maxfnS,4);
      XX                 = XJ;
      if ~SCALE
         XX(Reals,1)      = (XJ(Reals,1)-x_L(Reals))./x_D(Reals);
         XX(Reals,maxfnS) = (XJ(Reals,maxfnS)-x_L(Reals))./x_D(Reals);
      end
      FASTSTOP           = inf;
      SetvJ              = 0;
      fnStarJump         = inf;
      nTry               = 1;
      j                  = 1;
      while j < mTry-1
         j                 = j+1;
         fnStar            = fnS(j);
         if isnan(fJ(j,4))
            % No solution is computed for the current fnS(j)
            gnProb.CGO.fnStar = fnStar;
            % Use minimum on surface as initial x_0 for current target value
            gnProb.x_0        = min_sn_y;
            % Use previous solutions for all target values as matrix of initial points
            gnProb.X0         = XJ(:,[maxfnS,1:nTry]);
            gnProb.GO.ProbL.CGO = gnProb.CGO;
            gnR0              = gnSolve(gnProb);
            xTry              = gnR0.x_k(:,1);
            fTry              = gnR0.f_k;
            fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
            % Compute solution scaled to [0,1]
            XXTry            = xTry;
            if ~SCALE
               XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
            end
            % Accept new point and update
            nTry              = nTry+1;
            XJ(:,nTry)        = xTry;
            XX(:,nTry)        = XXTry;
            fJ(nTry,1:nfJ-DP) = fV;
            if DebugARBF
               fx             = fTrue(xTry,x_L,x_D,Prob,varargin{:});
               fJ(nTry,nfJ)   = fx;
            end
         else
            % A solution was previously computed for the current fnS(j)
            nTry              = nTry+1;
            fV                = fJ(nTry,1:nfJ-DP);
            xTry              = XJ(:,nTry);
            fTry              = fV(4);              % Same as fJ(nTry,4)
            XXTry             = XX(:,nTry);
         end
         doMNew            = fV(7);
         % Check if previous solutions was OK
         for k=nTry-1:2
            doMOld = fJ(k,7);
            if doMOld - doMNew > 1E-15 * doMNew
               if PriLev > 1
                  fprintf('    Must recompute solution %d\n',k);
               end
               % Use solution from new target value as initial point
               gnProb.x_0        = xTry;
               % Use best minimum in previous computation as additional initial value
               gnProb.X0         = XJ(:,k);
               % HKH0804 - must set correct fnStar
               gnProb.CGO.fnStar = fnS(k);
               gnProb.GO.ProbL.CGO = gnProb.CGO;
               gnR0              = gnSolve(gnProb);
               xTry              = gnR0.x_k(:,1);
               fTry              = gnR0.f_k;
               fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
               % Compute solution scaled to [0,1]
               XXTry            = xTry;
               if ~SCALE
                  XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
               end
               XJ(:,k)        = xTry;
               XX(:,k)        = XXTry;
               fJ(k,1:nfJ-DP) = fV;
               % Should also print out something
            end
         end
         D                 = normDist(XXTry,XX(:,nTry-1),dSQRT,IntVars,Reals);
         onB               = fV(2);
         if D(1) > 3*Gamma1 & onB <= minonB & onB == fJ(nTry-1,2)
            % Big difference between solutions, try to improve solution in case of bad optimum
            if PriLev > 1
               fprintf('    ');
               xprinti(fJ(1:nTry,2),'onB:');
               fprintf('    Recompute solution to possibly decrease D %f',D(1));
            end
            % Use solution from new target value as initial point
            gnProb.x_0        = xTry;
            % Use best minimum in previous computation as additional initial value
            gnProb.X0         = XJ(:,nTry-1);
            gnProb.CGO.fnStar = fnStar;
            gnProb.GO.ProbL.CGO = gnProb.CGO;
            gnR0              = gnSolve(gnProb);
            xTry              = gnR0.x_k(:,1);
            fTry              = gnR0.f_k;
            fV                = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals);
            % Compute solution scaled to [0,1]
            XXTry             = xTry;
            if ~SCALE
               XXTry(Reals) = (xTry(Reals)-x_L(Reals))./x_D(Reals);
            end
            XJ(:,nTry)        = xTry;
            XX(:,nTry)        = XXTry;
            fJ(nTry,1:nfJ-DP) = fV;
            D = normDist(XXTry,XX(:,nTry-1),dSQRT,IntVars,Reals);
            onB               = fV(2);
            if PriLev > 1
               % Should also print out something
               fprintf(' --- D now %f\n',D(1));
            end
         end
         doMOld            = fJ(nTry-1,7);
         fTryOld           = fJ(nTry-1,4);
         if PriLev > 1 & doMNew - doMInf > 1E-5 * doMInf
            fprintf('%2d: fnStar %15.7f doMInf %20.15f < doMNew %20.15f\n',j,fnStar,doMInf, doMNew);
            fprintf('  :               fTryInf %25.18f   fTryNew %25.18f\n',fTryInf,fTry);
         end
         if SetvJ & (doMNew - doMInf >= 1E-15 * doMInf | onB == d)
            if PriLev > 1
               fprintf('    ');
               fprintf('SetvJ: %2d: fNSTAR %15.7f doMInf %20.15f < doMNew %20.15f',j,fnStar,doMInf, doMNew);
               fprintf(' Decrease vJFac from %12.6f',vJFac);
            end
            vJVal = (min_sn - fnStar) / (vJFac*fRange);
            if RelLocErr < 1 & vJFac > 1E-3
               % if vJVal <= 1, vJFac =vJVal/vJ(end-2) * vJFac; end
               if vJVal <= 1, vJFac =max(vJFac/10,vJVal/vJ(end-2) * vJFac); end
               if PriLev > 1
                  fprintf(' to %12.6f',vJFac);
               end
            else
               if PriLev > 1
                  fprintf(' ! No Decrease, RelLocErr %f',RelLocErr);
               end
            end
            SetvJ = 0;
            if PriLev > 1
               fprintf('\n');
            end
         end
         if SetvJ & (nTry >= mTry-1)
            if vJFac < 1E3 & onB <= minonB
               vJFac = vJFac*10;
               if PriLev > 1
                  fprintf('%2d: Increase vJFac 10 times to %20.15f\n',j,vJFac);
               end
               SetvJ = 0;
             end
         end
         if PriLev > 1 & doMNew - doMInf > 1E-15 * doMInf
            fprintf('%2d: fNSTAR %15.7f doMInf %20.15f < doMNew %20.15f\n',j,fnStar,doMInf, doMNew);
         end
         %if doMOld - doMNew > 1E-15 * doMNew
         %   fprintf('%2d: fNSTAR %15.7f doMNew %20.15f < doMOld %20.15f\n',j,fnStar,doMNew, doMOld);
         %end
         if PriLev > 1 & doMOld - doMNew > 1E-5 * doMNew
            fprintf('%2d: fnStar %15.7f doMNew %20.15f < doMOld %20.15f\n',j,fnStar,doMNew, doMOld);
            fprintf('  :               fTryNew %25.18f   fTryOld %25.18f\n',fTry,fTryOld);
         end
         M = 0;
         if D(1) > 0 & mTry < maxfnS-1 & fnStar < fnStarJump
            if PriLev > 1
               fprintf('%2d: D %12.6f [%d]',j,D(1), onB)
               if D(1) > Gamma1 & onB == fJ(nTry-1,2), fprintf(' Too big range'); end
               fprintf('\n')
            end
            if D(1) > Gamma1 & onB <= minonB & onB == fJ(nTry-1,2)
               %HKH M = ceil(D(1)/Gamma1);
               M = 1+ceil(D(1)/Gamma1);
               if PriLev > 1
                  fprintf('    Must divide further, %d more points\n',M);
               end
               %HKH M = min([5,maxfnS-mTry-1,M]);
               M = min([7,maxfnS-mTry-1,M]);
               fnStarJump = fnStar;
            end
         end
         if 0 & onB > fJ(nTry-1,2) & D(1) > Gamma1 & fnStar < fnStarJump
            % HKH For LOCAL, skip making more fnStar if more components on bounds
            % New solution has more components on bounds
            if onB == fJ(nTry-1,2)+1 & D(1) > Gamma1 & onB < fJ(maxfnS,2)
               %HKH M = ceil(D(1)/Gamma1);
               M = 1+ceil(D(1)/Gamma1);
               if PriLev > 1
                  fprintf('    Must divide further II, %d more points\n',M);
               end
               %HKH M = min([5,maxfnS-mTry-1,M]);
               M = min([7,maxfnS-mTry-1,M]);
               fnStarJump = fnStar;
            elseif onB == fJ(nTry-1,2)+1
               % Normal change for new solution
            else
               % More than one component jumped to a bound
               if PriLev > 1
                  fprintf('    More than one component jumped to a bound: %d to %d \n',fJ(nTry-1,2),onB);
               end
               if nTry <= 10
                  M = min([10,maxfnS-mTry-1,3*(onB-fJ(nTry-1,2))]);
               else
                  M = min([5,maxfnS-mTry-1,onB-fJ(nTry-1,2)]);
               end
               fnStarJump = fnStar;
            end
         end
         if M > 0
            fnStarRange = fnS(nTry-1)-fnStar;
            %Newfn = fnS(nTry-1)-log([1:M])/log((M+1))*fnStarRange;
            Newfn = fnS(nTry-1)-[1:M]/(M+1)*fnStarRange;
            if PriLev > 1
               fprintf('    Set M = %d more points',M);
               fprintf(' fnStarRange %f',fnStarRange);
               fprintf(' fnStar %f',fnStar);
               fprintf(' fnSnTry %f',fnS(nTry-1));
               fprintf('\n');
               xprint(Newfn,'    Newfn:');
            end
            % Add new fnStar values first in list, move current to new position
            fnS  = [fnS(1:j-1),Newfn,fnS(j:end)];
            vJ   = [vJ(1:j-1),ones(1,M)*vJ(j),vJ(j:end)];
            % Move current fnStar solution to new position
            XJ(:,nTry+M)        = XJ(:,nTry);
            XX(:,nTry+M)        = XX(:,nTry);
            fJ(nTry+M,:)        = fJ(nTry,:);
            % Mark row nTry as not computed
            fJ(nTry,4)          = NaN;
            mTry                = mTry + M;
            nTry                = nTry - 1;
            j                   = j - 1;
         end
         % Avoid that range of fnStar values down to -inf with interesting solutions is not utilized
         if nTry == mTry-1 & onB ~= fJ(maxfnS,2)
            M                = maxfnS-mTry;
            mTry             = maxfnS;
            fnStarRange1     = fnS(1)-fnStar;
            fnStarRange2     = fnS(nTry-1)-fnStar;
            fnStarRange      = max(3*2^(-M)*fnStarRange1,fnStarRange2);
            Newfn            = fnS(nTry)-2.^[0:M-1]*fnStarRange;
            if PriLev > 1
               fprintf('    Increase Range to -inf:');
               fprintf(' fnStar %f',fnStar);
               fprintf(' fnStarRange1 %f',fnStarRange1);
               fprintf(' fnStarRange2 %f',fnStarRange2);
               fprintf(' fnStarRange %f',fnStarRange);
               fprintf('\n');
               xprint(Newfn,'    Newfn:');
            end
            % Add new fnStar values between current and -inf in list
            fnS  = [fnS(1:j),Newfn,fnS(end)];
            vJ   = [vJ(1:j),ones(1,M)*vJ(j),vJ(end)];
         end
         %if nTry == mTry-1 & onB == fJ(maxfnS,2) & onB == minonB
         if nTry == mTry-1 & onB == fJ(maxfnS,2)
            D  = normDist(XXTry,XX(:,maxfnS),dSQRT,IntVars,Reals);
            if PriLev > 1 & D(1) > Gamma1
               fprintf('    Here we should add some items down to -inf. D= %f',D(1));
               fprintf('\n');
            end
         end
         if FASTSTOP == nTry, break; end
         if FASTCODE & isinf(FASTSTOP) & ...
            (( onB > minonB | (abs(doMNew - doMInf) <= 1E-5 * doMInf)) & (onB <= fJ(maxfnS,2)))
            if M==0
               break;
            else
               FASTSTOP = nTry+M+1;
            end
         end
      end
      % Move -inf row to row nTry+1
      if nTry ~= maxfnS-1
         nTry              = nTry+1;
         XJ(:,nTry)        = XJ(:,maxfnS);
         XX(:,nTry)        = XX(:,maxfnS);
         fJ(nTry,:)        = fJ(maxfnS,:);
      end

      % Compute x_min solution scaled to [0,1]
      x_minS = x_min;
      if ~SCALE
         x_minS(Reals) = (x_minS(Reals)-x_L(Reals))./x_D(Reals);
      end
      D0                = normDist(x_minS,XX(:,1),dSQRT,IntVars,Reals);
      if PriLev > 1
         fprintf('Scaled distance from x_min to 1st solution (min_sn_y): D0 %f\n',D0(1));
      end

      % Print one row for each fnStar tried
      if PriLev > 2
         fnSPrint(XJ(:,1:nTry),fJ(1:nTry,:),min_sn,DebugARBF)
      end

      onBTry    = fJ(1:nTry,2);
      [g,DistGrp] = djCluster(XX(:,1:nTry),onBTry,epsX,PriLev > 2,IntVars,Reals);
      nGroup = g(end);

      fTry      = fJ(1:nTry,4);
      if DebugARBF
         fxTry      = fJ(1:nTry,end);
      else
         fxTry      = [];
      end
      doXTry    = fJ(1:nTry,6);
      % doMTry    = fJ(1:nTry,7);

      % xTry    = zeros(nGroup,1); % SEEMS NOT TO BE USED
      Gsize   = zeros(nGroup,1);
      GonBmin = zeros(nGroup,1);
      GonBmax = zeros(nGroup,1);
      ii      = 0;
      while ii < nGroup
         ii = ii+1;
         ix = find(g==ii);
         Gsize(ii) = length(ix);
         minonB      = min(onBTry(ix));
         GonBmin(ii) = minonB;
         GonBmax(ii) = max(onBTry(ix));
      end
      if PriLev > 1
         fprintf('RelGlobErr-1 %f ',RelGlobErr(1));
         fprintf('GlobErr-1 %f ',GlobErr(1));
         fprintf('\n');
         fprintf('RelGlobErr-2 %f ',RelGlobErr(2));
         fprintf('GlobErr-2 %f ',GlobErr(2));
         fprintf('\n');
         fprintf('RelLocErr    %f ',RelLocErr);
         fprintf('LocErr    %f ',LocErr);
         fprintf('\n');
         fprintf('\n');
      end
      xCand   = [];
      minonBG = min(GonBmin);
      % maxonBG = min(GonBmax); % NOT USED
      for i=1:nGroup
         j  = i;
         ix = find(g==j);
         iy = find(onBTry(ix) == GonBmin(j));
         iz = ix(iy);
         minDist  = min(DistGrp(iz));
         L  = max(DistGrp(iz))- minDist;
         % L  = max(DistGrp(ix));
         % l  = ceil(10*L); % NOTE l SEEMS NOT TO BE USED
         if j == 1
            M =  max(4,ceil(L/Gamma1));
            if D0(1) > L | D0(1) >= min(RelLocErr,Gamma1/3)
               D01=D0(1);
               k       = iz(1);
               xCand = [xCand,k];
               if PriLev > 1
                  fprintf('-------- WHY2? D01 %f ',D01)
                  fprintf('L %d ',L)
                  fprintf('M %d ',M)
                  fprintf('RelLocErr %e',RelLocErr)
                  fprintf('\n')
                  GrpPrint(i,j,k,'xCand1-DO',GonBmin,GonBmax,Gsize,fxTry,doXTry)
               end
            end
            % Generate M points, the points {ii/M*GrpDistance, i=1,M}
            for ii = 1:M-1
               %[vv,iy] = min(abs(DistGrp(iz)-(minDist+ii*L/M)));
               [vv,iy] = min(abs(DistGrp(iz)-(ii*L/M)));
               k       = iz(iy);
               if any(xCand==k) & length(iz) > iy
                  k       = iz(iy+1);
               end
               if ~any(xCand==k)
                  xCand = [xCand,k];
                  if PriLev > 1
                     GrpPrint(i,j,k,['xCand1-D',num2str(ii)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
                  end
               end
               %k = iz(ii);
               %if DistGrp(k)-D > 0.03
               %   xCand   = [xCand,k];
               %   D       = DistGrp(k);
               %   GrpPrint(i,j,k,'xCandmZ',GonBmin,GonBmax,Gsize,fxTry,doXTry)
               %end
            end
            k       = iz(end);
            xCand = [xCand,k];
            if PriLev > 1
               GrpPrint(i,j,k,['xCand1-D',num2str(M)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
            end
            %elseif j == nGroup
         elseif 0 & (j < nGroup & minonBG == GonBmin(i)) | (j == nGroup & minonBG == GonBmin(i))
            % Right now both inbetween groups 2,...,nGroup-1, as well as last group
            % elseif minonBG == GonBmin(i)
            M =  max(1,ceil(L/Gamma1));
            if M > 2
               k       = iz(1);
               xCand = [xCand,k];
               if PriLev > 1
                  GrpPrint(i,j,k,['xCand',num2str(i),'-1'],GonBmin,GonBmax,Gsize,fxTry,doXTry)
               end
            end
            % Generate M points, the points {ii/M*GrpDistance, i=1,M}
            for ii = 1:M-1
               %[vv,iy] = min(abs(DistGrp(iz)-(minDist+ii*L/M)));
               [vv,iy] = min(abs(DistGrp(iz)-(ii*L/M)));
               k       = iz(iy);
               if any(xCand==k) & length(iz) > iy
                  k       = iz(iy+1);
               end
               if ~any(xCand==k)
                  xCand = [xCand,k];
                  if PriLev > 1
                     GrpPrint(i,j,k,['xCand',num2str(i),'-',num2str(ii)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
                  end
               end
            end
            k       = iz(end);
            xCand = [xCand,k];
            if PriLev > 1
               GrpPrint(i,j,k,['xCand',num2str(i),'-',num2str(M)],GonBmin,GonBmax,Gsize,fxTry,doXTry)
            end
         end
      end
      ix = find(doXTry(xCand) > epsX);
      if isempty(ix)
         xCand = xCand(length(xCand));
         nCand = 1;
      elseif length(ix) ~= length(xCand)
         nCand  = length(xCand);
         if PriLev > 1
            fprintf('Remove %d candidates too close to X\n',nCand-length(ix))
         end
         xCand = xCand(ix);
         nCand  = length(xCand);
      else
         nCand  = length(xCand);
      end
      xLoc   = XJ(:,xCand);
      if nCand > 1
         if PriLev > 1
            xprinti(xCand,'xCand');
         end
         xC = xCand;
         for i = nCand:-1:2
            ii = xCand(i);
            for j=1:i-1
               jj = xCand(j);
               if tomsol(30,XX(:,jj),XX(:,ii))/dSQRT < 1E-6
                  xLoc = [xLoc(:,1:i-1),xLoc(:,i+1:end)];
                  xC   = [xC(1:i-1),xC(i+1:end)];
                  if PriLev > 1
                     fprintf('    Remove candidate %d, #%d, is too close to #%d\n',i,ii,jj);
                  end
                  break;
               end
            end
         end
         if PriLev > 1 & length(xC) < nCand
            xprinti(xC,'    Compressed xCand');
         end
         nCand    = size(xLoc,2);
         xCand    = xC;
         fLocSave = fTry(xCand(1:end));
         fnSSave  = fnS(xCand(1:end));
      end

      fLoc    = fTry(xCand(1));
      fnStar  = fnS(xCand(1));
      LocRed  = 0;   % Keep track if local grid points do reduce f(x)
      GlobRed = NaN; % No global grid now
      if nCand > 1
         if PriLev > 1
            fprintf('    Try %d new points. ',nCand);
         end
         modN = -4;
      elseif nCand == 1
         if PriLev > 1
            fprintf('    Try 1 new point.  ');
         end
         modN = -4;
      else
         if PriLev > 1
            fprintf('    No local found.  ');
         end
         fnStar = min_sn;
         modN   = -9;
      end
      if PriLev > 1
         fprintf('\n');
         if nCand > 0
            for i = 1:nCand
               %fprintf('xNew #%d:',xCand(i));
               %fprintf('\n');
               xprint(xLoc(:,i),'  xNew:');
            end
         end
      end
      %keyboard

   else
      error('Illegal modN')
   end



   % No point in not accepting local solution if Pure Integer problem
   if length(IntVars) == d, ucOK = 1; end

   % ********** MINIMIZE gn(y) **********
   if modN >= 0 & Iter > 1
      [v,ix] = find(Its.modN(1:Iter-1)==modN);
      if ~isempty(ix)
         ixlast = max(ix);
         if Its.FLOWER(ixlast) == Its.nFunc(ixlast)
            PROGRESS(1+modN) = 1;
         else
            PROGRESS(1+modN) = 0;
         end
         if PriLev > 1 & modN == N
            fprintf('Last local step %d. ',ixlast);
            if PROGRESS(N+1)
               fprintf('FUNCTION REDUCTION. ');
            else
               fprintf('NO FUNCTION REDUCTION. ');
            end
            fprintf('RelLocErr    %f ',RelLocErr);
            fprintf('LocErr    %f',LocErr);
            fprintf('\n');
         end
      end
   end
   if (modN == N) & ucOK
      if PriLev > 2
         fprintf('Local search OK, not minimizing gn_f\n')
      end
      if OLDuc == 0
         xGlob = [];
         fGlob = [];
         xLoc  = min_sn_y;
         fLoc  = min_sn;
      elseif OLDuc == 1
         % Minimize s_n(y) using local optimization
         %  -  Solve  min[s_n(y)] with a local optimizer by starting from
         %     the interpolation point with the least function value
         % i.e. Find a local solution close to (x_min,fMin)
         % Set parameters used in snSolve
         snProb.CGO.globalSolver = GOlocalSolver;
         snProb.CGO.localSolver  = GOlocalSolver;
         snProb.CGO.fnStar       = NaN;
         snProb.x_0              = x_min;
         %HKH NEW
         snProb.PriLev           = PriSub;
         snProb.CGO.modN         = modN;
         snProb.GO.ProbL.CGO     = snProb.CGO;
         snResult                = snSolve(snProb);
         fLoc                    = snResult.f_k;
         xLoc                    = snResult.x_k(:,1);
         if PriLev > 1
            xprint(fLoc,'fLoc  ');
            xprint(min_sn,'min_sn');
            xprint(xLoc,'xLoc    ');
            xprint(min_sn_y,'min_sn_y');
         end
         Dist1 = tomsol(30,xLoc,x_min);
         Dist2 = tomsol(30,min_sn_y,x_min);
         if Dist1 > Dist2
            if PriLev > 1
               fprintf('USE min_sn_y instead!!! ')
               fprintf('Dist1 %f ', Dist1);
               fprintf('Dist2 %f ', Dist2);
               fprintf('\n');
            end
            xLoc  = min_sn_y;
            fLoc  = min_sn;
         end
      else
         % L-step Case 1
         snProbD1.user.X      = X;
         snProbD1.c_L(dCon(1)+1) = (5E-4*max(1E-3,norm(x_min)))^2;
         snProbD1.ixDist      = fIdx;
         snProbD1.x_0         = x_min;
         snProbD1.mNonLin     = dCon(1)+1;
         snProbD1.CGO.fnStar  = NaN;
         %DbetaV               = [0.005,0.01:0.01:0.05,0.10:0.05:0.50];
         DbetaV               = [0.005, 0.010, 0.012, 0.015];
         xLoc                 = [];
         fLoc                 = [];
         if DebugARBF
            if SCALE
               O_new   = tomsol(9, x_L, min_sn_y, x_D);
            else
               O_new   = min_sn_y;
            end
            [f00OrgLoc,CcNew] = cgo_fc(O_new, Prob, varargin{:});
            fOrgLoc           = cgo_pf(f00OrgLoc, CcNew, Prob);
            %NLP_x = [];
            %fOrgLoc  = nlp_f(O_new, Prob, varargin{:});
            if PriLev > 1
               fprintf('   RBF min_sn               ');
               fprintf('fx %11.8f ',fOrgLoc);
               fprintf('f %f ', min_sn);
               fprintf('            ');
            end
            [onB_sn, doX_sn, doM_sn] = statGN(snProb.x_L,snProb.x_U,min_sn_y,...
               x_min,X,epsX,IntVars,Reals);
            if PriLev > 1
               fprintf('[%d] ',onB_sn);
               if isempty(IntVars)
                  fprintf('doX %6.3f ',doX_sn);
                  fprintf('doM %8.5f ',doM_sn);
                  fprintf('            ');
               else
                  fprintf('doX %6.3f/%d ',doX_sn);
                  fprintf('doM %8.5f/%d ',doM_sn);
                  fprintf('            ');
               end
            end
            if ~isempty(xOptS)
               doO_sn = min(tomsol(30,min_sn_y,xOptS));
            else
               doO_sn = [];
            end
            if PriLev > 1
               fprintf('doO %8.5f ',doO_sn);
               xprint(min_sn_y,'x:',' %12.8f',10)
            end
         end
         for i = 1:length(DbetaV)
            Dbeta = DbetaV(i);
            snProbD1.c_U(dCon(1)+1) = (Dbeta*pDist)^2;
            if strcmpi(localSolver,globalSolver)
               snR  = tomRun(GOlocalSolver,snProbD1,PriLev-4);
            else
               snR  = tomRun(localSolver,snProbD1,PriLev-4);
            end
            x_k             = snR.x_k;
            [onB, doX, doM] = statGN(snProbD1.x_L,snProbD1.x_U,x_k,x_min,X,epsX,IntVars,Reals);
            doS             = min(tomsol(30,x_k,min_sn_y));
            if ~isempty(xOptS)
               doO = min(tomsol(30,x_k,xOptS));
            else
               doO = [];
            end
            %if i==6 % Save the solution at 0.05
            %   xLoc     = x_k;
            %   fLoc     = snR.f_k;
            %   snResult = snR;
            %end
            %if i == length(DbetaV)
            %   DbetaOpt = DbetaV(6);
            %else
            %   DbetaOpt = Dbeta;
            %end
            c_k = snR.c_k(dCon(1)+1);
            if isempty(snR.v_k)
               % [v_k,Zv,P,Z_L,cErr,ceq,cineq,gProj] = LagMult(snProbD1,snR);
               % [v_k] = LagMult(snProbD1,snR);
               LagMul   = c_k >= snProbD1.c_U(end)-cTol;
            else
               LagMul   = snR.v_k(end);
            end
            if SCALE
               O_new   = tomsol(9, x_L, x_k, x_D);
            else
               O_new   = x_k;
            end
            if DebugARBF
               [f00OrgLoc,CcNew] = cgo_fc(O_new, Prob, varargin{:});
               fOrgLoc           = cgo_pf(f00OrgLoc, CcNew, Prob);
               %NLP_x = [];
               %fOrgLoc  = nlp_f(O_new, Prob, varargin{:});
            end
            if LagMul == 0 & min(doX,1) > 1E-7 * norm(x_min)
               xLoc     = x_k;
               fLoc     = snR.f_k;
               doXS     = doX;
               doMS     = doM;
               doSS     = doS;
               if PriLev > 1
                  fprintf('OK Db %f c_k %f ', Dbeta, c_k);
                  if DebugARBF
                     fprintf('fx %11.8f ',fOrgLoc);
                  end
                  fprintf('f %f ', snR.f_k);
                  fprintf('LagM %6.1f ',LagMul);
                  %fprintf('ucOK %d ',ucOK);
                  fprintf('[%d] ',onB);
                  if isempty(IntVars)
                     fprintf('doX %6.3f ',doX);
                     fprintf('doM %8.5f ',doM);
                     fprintf('doS %7.5f ',doS);
                  else
                     fprintf('doX %6.3f/%d ',doX);
                     fprintf('doM %8.5f/%d ',doM);
                     fprintf('doS %7.5f/%d ',doS);
                  end
                  if ~isempty(xOptS)
                     fprintf('doO %8.5f ',doO);
                  end
                  xprint(xLoc,'x:',' %12.8f',10)
                  %fprintf('\n');
               end
               snResult = snR;
               break;
            else
               if PriLev > 1
                  %PrintResult(snR,2);
                  fprintf('   Db %f c_k %f ', Dbeta, c_k);
                  if DebugARBF
                     fprintf('fx %11.8f ',fOrgLoc);
                  end
                  fprintf('f %f ', snR.f_k);
                  fprintf('LagM %6.1f ',LagMul);
                  %fprintf('ucOK %d ',ucOK);
                  fprintf('[%d] ',onB);
                  if isempty(IntVars)
                     fprintf('doX %6.3f ',doX);
                     fprintf('doM %8.5f ',doM);
                     fprintf('doS %7.5f ',doS);
                  else
                     fprintf('doX %6.3f/%d ',doX);
                     fprintf('doM %8.5f/%d ',doM);
                     fprintf('doS %7.5f/%d ',doS);
                  end
                  if ~isempty(xOptS)
                     fprintf('doO %8.5f ',doO);
                  end
                  xprint(x_k,'x:',' %12.8f',10)
                  %fprintf('\n');
               end
            end
         end
         snProbD1.PriLevOpt = 0;
         snProbD1.xInit     = 20;
         snProbD1.xEqTol    = 1E-4;
         %DbetaV            = [0.005,0.010];
         DbetaV             = [];
         for i = 1:length(DbetaV)
            Dbeta = DbetaV(i);
            snProbD1.c_U(dCon(1)+1) = (Dbeta*pDist)^2;
            snR  = tomRun('multiMin',snProbD1,PriLev-4);
            x_k             = snR.x_k(:,1);  % HKH fix: For now set 1st x solution
            [onB, doX, doM] = statGN(snProbD1.x_L,snProbD1.x_U,x_k,x_min,X,epsX,IntVars,Reals);
            doS             = min(tomsol(30,x_k,min_sn_y));
            if ~isempty(xOptS)
               doO = min(tomsol(30,x_k,xOptS));
            else
               doO = [];
            end
            %if i==6 % Save the solution at 0.05
            %   xLoc     = x_k;
            %   fLoc     = snR.f_k;
            %   snResult = snR;
            %end
            %if i == length(DbetaV)
            %   DbetaOpt = DbetaV(6);
            %else
            %   DbetaOpt = Dbeta;
            %end
            c_k = snR.c_k(dCon(1)+1);
            if isempty(snR.v_k)
               % [v_k,Zv,P,Z_L,cErr,ceq,cineq,gProj] = LagMult(snProbD1,snR);
               % [v_k] = LagMult(snProbD1,snR);
               LagMul   = c_k >= snProbD1.c_U(end)-cTol;
            else
               LagMul   = snR.v_k(end);
            end
            if SCALE
               O_new   = tomsol(9, x_L, x_k, x_D);
            else
               O_new   = x_k;
            end
            if DebugARBF
               [f00OrgLoc,CcNew] = cgo_fc(O_new, Prob, varargin{:});
               fOrgLoc           = cgo_pf(f00OrgLoc, CcNew, Prob);
               %NLP_x = [];
               %fOrgLoc  = nlp_f(O_new, Prob, varargin{:});
            end
            if PriLev > 1
               fprintf('   Db %f c_k %f ', Dbeta, c_k);
               if DebugARBF
                  fprintf('fx %11.8f ',fOrgLoc);
               end
               fprintf('f %f ', snR.f_k);
               fprintf('LagM %6.1f ',LagMul);
               fprintf('[%d] ',onB);
               if isempty(IntVars)
                  fprintf('doX %6.3f ',doX);
                  fprintf('doM %8.5f ',doM);
                  fprintf('doS %7.5f ',doS);
               else
                  fprintf('doX %6.3f/%d ',doX);
                  fprintf('doM %8.5f/%d ',doM);
                  fprintf('doS %7.5f/%d ',doS);
               end
               if ~isempty(xOptS)
                  fprintf('doO %8.5f ',doO);
               end
               xprint(x_k,'x:',' %12.8f',10)
            end
         end
         snProbD1.PriLevOpt = DebugPriLev;
      end
      if isempty(xLoc)
         doXS     = doX;
         doMS     = doM;
         doSS     = doS;
         xGlob = [];
         fGlob = [];
         if RelLocErr < 0.01 | doM_sn(1) > doMS(1)
            if PriLev > 1
               fprintf('PICK MIN ON RBF SURFACE %15.8f\n',min_sn);
            end
            xLoc  = min_sn_y;
            fLoc  = min_sn;
         else
            snResult = snR;
            xLoc  = snResult.x_k;
            fLoc  = snResult.f_k;
            if PriLev > 1
               fprintf('TAKE LOCAL SOLUTION AS FAR AS POSSIBLE FROM fMin, doM %8.5f\n',doMS);
            end
         end
      else
         if (RelLocErr <  0.01 & doM_sn(1) <  doMS(1)) | ...
               (RelLocErr >= 0.01 & doM_sn(1) >= doMS(1))
            if PriLev > 1
               if RelLocErr < 0.01
                  fprintf('LOCAL PHASE,  RelLocError %10.5f. ',RelLocErr);
               else
                  fprintf('GLOBAL PHASE, RelLocError %10.5f. ',RelLocErr);
               end
               fprintf('USE MIN ON RBF SURFACE %15.8f\n',min_sn);
            end
            xLoc  = min_sn_y;
            fLoc  = min_sn;
         end
      end
   end
   if modN == -2
      % Use min_sn, minimum on RBF surface
      xLoc     = snResult.x_k;
      fLoc     = snResult.f_k;
   end
   if modN == -8 | modN == -9 | modN == -10
      % Use min_sn, minimum on RBF surface
      xLoc     = snResult.x_k;
      fLoc     = snResult.f_k;
   end
   if modN == -5
      % Use min_sn, minimum on RBF surface
      xLoc     = snResult.x_k;
      fLoc     = snResult.f_k;
   end
   if modN == -6
      % Use min_sn, minimum on RBF surface and a NEWTON step
      % global NLP_x NLP_f NARG
      NLP_x=[]; NLP_f=[]; NARG = [];
      sn_fx = nlp_f(x_min,snProb);
      sn_gx = nlp_g(x_min,snProb);
      sn_Hx = nlp_H(x_min,snProb);
      p    = -sn_Hx\sn_gx;
      if PriLev > 1
         fprintf('sn_f %11.6f ',sn_fx);
         fprintf('sn_f %11.6f (equal???) ',snNew);
         xprint(p,'p');
      end
      xLoc = max(x_LL,min(x_UU,snResult.x_k+p));
      fLoc = snResult.f_k;
   end
   if modN == -7
      % Use min_sn, minimum on RBF surface
      xLoc     = snResult.x_k(:,1); %HKH Fix, for now use 1st global opt point
      fLoc     = snResult.f_k;
   end
   %if modN ~= N | (modN == N & ~ucOK)
   if (modN == N & ~ucOK) | modN == -1 | modN == -11
      if modN == -1 | modN == -11
         gnProb.CGO.fnStar = 0;
         gnProb.CGO.alpha  = 0;
      end
      gnProb.CGO.modN = modN;
      if PriLev > 1 & fnStar > fMin
         fprintf('------ WARNING!!! fnStar > fMin. ');
         fprintf('fnStar %18.12f ',fnStar);
         fprintf('fMin %18.12f ',fMin);
         fprintf('min_sn %18.12f ',min_sn);
         fprintf('\n');
         fnStar = fMin - 0.01*abs(fMin);
      end

      % Set iteration dependent parameters used in gnSolve, output from snSolve
      gnProb.IX   = snResult.multiMin.IX;
      gnProb.xOpt = snResult.multiMin.xOpt;

      gnProb.GO.ProbL.CGO = gnProb.CGO;
      gnResult = gnSolve(gnProb);
      xLoc     = gnResult.xLoc(:,1); %HKH Fix, for now use 1st global opt point
      fLoc     = gnResult.fLoc;
   end

   % Best point found on surface is xLoc - set as xNew
   if modN == -3 | modN == -4
      if size(xLoc,2) == 1
         % Last point in set
         xSave    = [];
      else
         fLocSave = fLocSave(2:end);
         fnSSave  = fnSSave(2:end);
         xSave    = xLoc(:,2:end);
         xLoc     = xLoc(:,1);
      end
   end
   if size(xLoc,2) > 1
      'NOT POSSIBLE TO REACH THIS LINE - BIG ERROR - MORE THAN 1 OPTIMAL SOLUTION'
      modN
      % More than one new point selected
      modN  = -3-modN;
      xSave = xLoc(:,2:end);
      xLoc  = xLoc(:,1);
   end
   xNew = xLoc;
   ix   = find(all(X==xNew*ones(1,size(X,2))));
   if ~isempty(ix)
      % New point is already sampled in a previous step
      RESCUE = 1;
      if PriLev > 0
         'RESCUE'
      end
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
   gnProb.CGO.alpha  = 0;

   % *************** UPDATE ***************
   % Remove information for response surface problem
   % global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];

   % Distance between new point and best point found or closest point in X
   [onB, doX, doM] = statGN(gnProb.x_L,gnProb.x_U,xNew,x_min,X,epsX,IntVars,Reals);
   % Distance between surface minimum and best point found or closest point in X
   [onB_sn, doX_sn, doM_sn] = statGN(snProb.x_L,snProb.x_U,min_sn_y,...
      x_min,X,epsX,IntVars,Reals);
   % Distance between surface minimum and known optimal point, if given
   if ~isempty(xOptS)
      doO_sn = min(tomsol(30,min_sn_y,xOptS));
   else
      doO_sn = [];
   end
   % Distance between new point and known optimal point, if given
   if ~isempty(xOptS)
      doO = min(tomsol(30,xNew,xOptS));
   else
      doO = [];
   end
   % Distance between new point and min on RBF surface
   doS = min(tomsol(30,xNew,min_sn_y));

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
         [f00New,CcNew] = cgo_fc(O_new, Prob, varargin{:});
         fNew           = cgo_pf(f00New, CcNew, Prob);
         % fNew  = nlp_f(O_new, Prob, varargin{:});
         nFunc = nFunc + 1; % Total number of sampled points
         fPen  = fNew;
         if dLin > 0
            L = Prob.A*O_new;
            fPen = fPen+sum(max(0,max(Prob.b_L-bTol*max(1,abs(Prob.b_L))-L,...
                                      L-bTol*max(1,abs(Prob.b_U))-Prob.b_U)));
         end
         if dCon(1) > 0 
            C    = cgo_c(O_new, Prob, varargin{:});
            fPen = fPen+sum(max(0,max(Prob.c_L-cTol*max(1,abs(Prob.c_L))-C,...
	                              C-cTol*max(1,abs(Prob.c_U))-Prob.c_U)));
            nCon = nCon + 1;
         end
      else
         Update = -4;
         fNew   = NaN;
      end
   end
   time  = fix(clock);

   if ~isnan(fNew)
      Dx     = tomsol(30,X,xNew);
      L      = abs(F(1:n-1)-F(n))./Dx(1:n-1,1);
      LipLow = min(LipLow,min(L));
      LipUpp = max(LipUpp,max(L));
   end

   if(usecgolib == 1 & Update == 1)
     % Add the point. Need checks on:
     %  * Close point
     %  * Identical point
     cgolib(102, Prob.CGOLIB.obs, xNew, Prob.CGOLIB.datacol, fNew(1));
     % Recompute the transforms. This should be improved, so
     % that transforms that could be updated row-wise will
     % only update the changed row and not all rows.
     cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
     % The surface needs to be interpolated again.
     % Need checks on:
     %  * Interplation ok?
     %  * And if we have update, we need to have mechanism
     %    to recompute LU and so on if interpolation was
     %    not ok.
     cgolib(302, Prob.CGOLIB.rbf);
     % For now, assume everything was ok
     Update = 1;
     
     % F_m0 should be set if REPLACE == 1. 
     if(REPLACE >= 1)
       nrows = cgolib(107, Prob.CGOLIB.obs);
       F_m0 = cgolib(110, Prob.CGOLIB.obs, nrows, 1, [], Prob.CGOLIB.rbfcol);
     end
     F0   = [F;fNew];
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
          if PriLev > 1
             fprintf('Restriced f values, indices:')
             xprinti(ix);
          end
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
      Fpen     = [Fpen;fPen];
      F00      = [F00;f00New];
      if ~isempty(CcNew)
         Cc    = [Cc,CcNew];
      end
      X        = [X xNew];
      O        = [O O_new];
      O_pre    = O_new;
      n        = n + 1;
      % Is the function value in the new point less than fMin?
      % If feasible, compare fNew,fMin, otherwise fPen and fMin
      % Or now feasible, but previously not?
      fRed  = fMin - fPen;
      if (Feasible & fNew < fMin & fPen==fNew) | ...
            (~Feasible & fPen < fMin) | (~Feasible & fPen == fNew)
         Feasible = fNew == fPen;
         fMin     = fPen;
         CcMin    = CcNew;
         fMinIter = Iter;
         fIdx     = length(Fpen);
         x_min    = xNew;
         O_min    = O_new;
         FLOWER   = nFunc;
         % FLOW     = nFunc; % NOT USED NOW
      end
      NOUPDATE = 0;
   elseif Update == -2
      % Update == -2 New point bad, even refactorization failed
      % Infeasible problem, ill-conditioning?

      Dist = tomsol(30,X(:,fIdx),X);
      Dist(fIdx) = Inf;

      control = -1;
      VALUE = 0;
      while control < 0
	if(usecgolib == 1)
	  error('This can''t happen for cgolib at the moment.');
	else
	  tomsol(25) % Deallocates memory
	  [minDist ixDist] = min(Dist);
          if PriLev > 1
	     fprintf('Minimal distance to X set');
	     fprintf(' %s',minDist);
	     fprintf('\n');
          end
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
	  F00  = F00(ix);
          if ~isempty(CcNew)
             Cc    = Cc(:,ix);
          end
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
	  disp('make init again')
	  control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
	  VALUE = 0;
	end
      end
      if control < 0
         fprintf('New initial interpolation failed');
         Result.ExitFlag  = 1;
         Result.ExitText  = 'New Initial interpolation failed';
         Result.CGO.fGoal = fGoal;
         Result.CGO.Its   = Its;
         Result.DIGIT     = TESTDIG;
         Result           = endSolve(Prob,Result);
         tomsol(25) % Deallocates memory
         return     %Something is really wrong
      end
      NOUPDATE = 0;
   elseif Update == -3 | Update == -4 | Update == -5
      % Update == -3 Point too close to old point
      % Update == -4 Point identical to old point
      % Update == -5 No feasible point found on surface
      NOUPDATE = NOUPDATE+1;
      SUCCESS   = 0;
      if PriLev > 1
         fprintf('!!!!! Update impossible\n')
      end
   end
   if RESCUE > 0
      SUCCESS   = 0;
   end
%HKHDEBUGprintmat(X')

   if isempty(fGoal) | isinf(fGoal)
      RelErr = NaN;
   else
      if fGoal == 0
         RelErr = abs(fMin-fGoal);
      else
         RelErr = abs(fMin-fGoal)/abs(fGoal);
      end
   end
   surfErr = fNew-snNew;
   if modN == N | modN == -8 | modN == -9 | modN == -10 | modN == -1 | ...
         modN == -5 | modN == -6 | modN == -7 | modN == -11
      LocErr = surfErr;
      if fNew ~= 0
         RelLocErr = abs(LocErr/fNew);
      else
         RelLocErr = abs(LocErr);
      end
   end
   if modN == 0
      j = modN + 1;
      GlobErr(j) = surfErr;
      if fNew ~= 0
         RelGlobErr(j) = abs(GlobErr(j)/fNew);
      else
         RelGlobErr(j) = abs(GlobErr(j));
      end
      if FLOWER  == nFunc;
         % Check if we improved local convergence
         if fNew ~= 0
            if abs((Its.fMin(Iter-1)-fNew)/fNew) < 0.01;
               LocErr = surfErr;
               RelLocErr = abs(LocErr/fNew);
            end
         else
            if abs(Its.fMin(Iter-1)-fNew) < 0.01;
               LocErr = surfErr;
               RelLocErr = abs(LocErr);
            end
         end
      end
   end
   if modN == -3
      if fRed > 1E-6*max(1,abs(fMin))
         GlobRed = 1;
      end
   end
   if modN == -4
      if fRed > 1E-6*max(1,abs(fMin))
         LocRed = 1;
      end
   end

   Its.Iter(Iter)             = Iter;
   Its.n(Iter)                = n;
   Its.nFunc(Iter)            = nFunc;
   Its.modN(Iter)             = modN;
   Its.fnStar(Iter)           = fnStar;
   Its.fMin(Iter)             = fMin;
   Its.fRed(Iter)             = fRed;
   Its.RelErr(Iter)           = RelErr;
   Its.surfErr(Iter)          = surfErr;
   Its.FLOWER(Iter)           = FLOWER;
   Its.fMinIter(Iter)         = fMinIter;
   Its.fNew(Iter)             = fNew;
   Its.minSn(Iter)            = min_sn;
   Its.onBminSn(Iter)         = onB_sn;
   Its.distminSn2X(Iter,:)    = doX_sn;
   Its.distminSn2xMin(Iter,:) = doM_sn;
   Its.onBxNew(Iter)          = onB;
   Its.distxNew2X(Iter,:)     = doX;
   Its.distxNew2xMin(Iter,:)  = doM;
   Its.distxNew2minSn(Iter,:) = doS;
   if ~isempty(xOptS)
      % global NLP_x NLP_f NARG
      NLP_x=[]; NLP_f=[]; NARG = [];
      snOptf             = nlp_f(xOptS,snProb);
      snOptg             = norm(nlp_g(xOptS,snProb));
      if isnan(snOptg)
         snOptH = [];
         snOptE = NaN;
      else
         snOptH = nlp_H(xOptS,snProb);
         snOptE = sum(eig(snOptH)<0);
      end
      dXO                = min(tomsol(30,xOptS,X));

      Its.snOptf(Iter)   = snOptf;
      Its.snOptg(Iter)   = snOptg;
      Its.snOptE(Iter)   = snOptE;
      Its.dXO(Iter)      = dXO;
   end

   % New variables saved
   Its.LipUpp(Iter)      = LipUpp;
   Its.LipLow(Iter)      = LipLow;
   Its.minSn_y(:,Iter)   = min_sn_y;
   Its.xInf(:,Iter)      = xInf;
   Its.fInf(Iter)        = fInf;
   Its.snNew(Iter)       = snNew;

   if SAVE == 1
      fMinIdx  = fIdx(1);
      rngState = rand('state');
      saveCGO(1,'arbfMIP',Name,O,F,X,F_m,F00,Cc,nInit,Fpen,fMinIdx,rngState);
   end
   

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
      fprintf(' fnStar%7.3f',fnStar);
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
      fprintf(' snErr %7.4f', surfErr);
      fprintf('\n');
      if ~isempty(CcMin)
         fprintf('CcMin:  ');
         fprintf('%11.8f ', CcMin);
      end
      xprint(O_new,'  xNew:',' %12.8f',8)
      %HKH Used in rbfSolve ???
      %fprintf(' fLoc %11.8f', fLoc);
      %NEWHKHfprintf(' %11.8f', fGlob);
      %fprintf(' min_sn %7.4f', min_sn);

      if PriLev > 2
         fprintf('  min_sn: %7.4f ',min_sn);
         fprintf('[%d] ',onB_sn);
         if isempty(IntVars)
            fprintf('doX %7.4f ',doX_sn);
            fprintf('doM %7.4f ',doM_sn);
         else
            fprintf('doX %f7.4/%d ',doX_sn);
            fprintf('doM %f7.4/%d ',doM_sn);
         end
         if ~isempty(xOptS)
            fprintf('doO %7.4f ',doO_sn);
         end
         fprintf('xNew: [%d] ',onB);
         if isempty(IntVars)
            fprintf('doX %7.4f ',doX);
            fprintf('doM %7.4f ',doM);
         else
            fprintf('doX %7.4f/%d ',doX);
            fprintf('doM %7.4f/%d ',doM);
         end
         if ~isempty(xOptS)
            fprintf('doO %7.4f ',doO);
         end
         fprintf('doS %7.4f ',doS);
         fprintf('\n');
      end

      if PriLev > 3
         fprintf('  snNew-min_sn %f ',snNew-min_sn);
         if idea == 1
            fprintf('snNew-fnStar %f ',snNew-fnStar);
         end
         fprintf('snNew-fNew %f ',snNew-fNew);
         if abs(myNew) > 1E8
            fprintf('myNew %e ',myNew);
         else
            fprintf('myNew %f ',myNew);
         end
         %HKH Used in rbfSolve
         %fprintf('fRed %f. ',fRed)
         %if modN ~= N % | (Rescue == 0 & ~ucOK)
         %   %fprintf('hn %f ',-1/gnResult.f_k);
         %   fprintf('hn %f ',-1/gnResult.f_k);
         %   if idea == 1
         %      %fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+1/gnResult.f_k);
         %      fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+1/gnResult.f_k);
         %   end
         %end
         fprintf('\n');
      end
      if ~isempty(xOptS) & PriLev > 3
         fprintf('  dXO (min||xOpt-X||) %f ',dXO);
         fprintf('sn_f(xOpt) %8.5f ',snOptf);
         fprintf('||sn_g(xOpt)|| %f ',snOptg);
         fprintf('negeig(sn_H(xOpt)) %d ',snOptE);
         fprintf('\n');
      end
      if PriLev > 2
         fprintf('  GlobRed %d. ',GlobRed)
         fprintf('LocRed %d. ',LocRed)
         fprintf('fRed %f. ',fRed)
         fprintf('LipUpp %f.',LipUpp)
         fprintf('LipLow %f. ',LipLow)
         fprintf('\n');
      end

   end

   % ********** CONVERGENCE TEST **********
   if n >= nMax, convflag = 7; end
   if convflag == 0
      convflag = CGOisClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
   end
   if convflag == 0
      if NOUPDATE > MaxCycle
         convflag = 4;
      elseif SAME1 > MaxCycle
         convflag = 5;
      elseif SAME2 > MaxCycle
         convflag = 6;
      end
   end

   % ------------------------------------------
   % Determine what type of step next iteration
   % ------------------------------------------
   if modN == -8
      if GlobRed | SUCCESS == 0
         % Redo global step if successful global search, after min_sn used
         modN = -1;
      else
         % Turn to local search if NOT successful global search
         modN = 0;
      end
      distxx =norm(x_min-min_sn_y)/max(1,norm(x_min));
      if PriLev > 1 & ...
       (RelLocErr/max(1E-300,distxx) < 0.0001 | (RelLocErr < 1E-4 & distxx==0)) 
            fprintf('NOTE - CONVERGENCE TO LOCAL MINIMUM in Global Grid, %17.10e\n',...
              RelLocErr/distxx);
      end
   elseif modN == -9
      if LocRed & SUCCESS
         % Redo local step if successful local search, after min_sn used
         % HKH ???
         modN = 0;
         %modN = -1;
      else
         % Turn to global search if NOT successful local search
         modN = -1;
      end
      if RelLocErr < 0.0001
         % Take a standard local step next step
         %modN = N-1;
         %modN = -6;
         %modN = N-1;
         modN = 10;
      end
      distxx =norm(x_min-min_sn_y)/max(1,norm(x_min));
      if PriLev > 1 & ...
       (RelLocErr/max(1E-300,distxx) < 0.0001 | (RelLocErr < 1E-4 & distxx==0))
         fprintf('NOTE - CONVERGENCE TO LOCAL MINIMUM in Local Grid, %17.10e\n',...
              RelLocErr/distxx);
      end
   elseif modN == -10
      % Turn to global search always
      modN = -1;
   elseif modN == -5
      if LocRed
         % Redo local step
         modN = -5;
      else
         % Turn to local search if NOT successful local step
         modN = 0;
      end
   elseif modN == -6
      if fRed > 0
         % Newton step was successful
         modN = -6;
      else
         modN = -7;
      end
   elseif modN == -7
      if LocRed
         % Go back to Newton step
         modN = -6;
      else
         % Turn to local search if NOT successful local step
         modN = 0;
      end
   elseif modN == -3 & isempty(xSave) & AddinfStep
      % Last point in set, now add infStep
      % modN     = abs(modN)-3;
      modN       = -11;
      AddinfStep = 0;
   elseif (modN == -3 & isempty(xSave)) | modN == -11
      % Last point in set, after possibly adding infStep
      % modN     = abs(modN)-3;
      modN = -8;
   elseif modN == -4 & isempty(xSave)
      % Last point in set
      modN = -9;
   end


   % ------------------------------------------
   % Print digits of convergence to given fGoal
   % ------------------------------------------
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

   if ExDTEST
      %NHQ ExD test
      % Distance to x_opt from current fMin.
      minXdist(:,Iter) = [min( tomsol(30,x_opt,O(:,fIdx)) ) ; fMin ; O(:,fIdx)];

      % X_1% norm, any point inside box centered around x_opt?
      while any( all(x_opt_L <= repmat(O(:,end),1,nrXopt)) & ...
            all(x_opt_U >= repmat(O(:,end),1,nrXopt))  )
         %all(x_opt_L <= O(:,end))  &&  all(x_opt_U >= O(:,end))
         xDIGITS = [xDIGITS [Iter nFunc]'];
         %x_box   = 0.1*x_box;
         x_box    = x_box/sqrt(10);
         x_opt_L = x_opt - repmat(0.5*x_box,1,nrXopt);
         x_opt_U = x_opt + repmat(0.5*x_box,1,nrXopt);
      end
   end

   % -------------- Result saving -------------------
end

% *******************************************************
% *************** END MAIN ITERATION LOOP ***************
% *******************************************************

%NHQ
if RMSError & d <= 4 
   % Calculate RMS error of the interpolated surface.
   % Performed for d = 2,3 or 4.
   % Crossvalidate surface before setting free cgolib.
   fprintf(' Starting valuation of the final RBF-surface.\n')
   % Always use original values in order to compare with real f(x).
   DACEProb.CGOLIB.TRANSFORM = 0;
   DACEProb.GO.ProbL.CGOLIB.TRANSFORM = 0;

   % Number of sampled points in RMS-evaluation
   nS = 0;
   if d == 2
      xGrid = linspace(0,1,41);
      N = length(xGrid)^d;
      sFx = zeros(N,1);
      fI  = true(N,1);        %Feasibility Index
      kk = 0;
      for ii = xGrid
         for jj = xGrid
            kk = kk+1;
            x  = [ii jj]';
            xF = x_L + x.*x_D;
            if dLin + dCon(1) > 0
               cVal = 0;
               if dLin > 0
                  L = Prob.A*xF;
                  cVal = cVal+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
               end
               if dCon(1) > 0
                  C = cgo_c(xF, Prob, varargin{:});
                  nCon = nCon + 1;
                  cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
               end
               if cVal > 0
                  fI(kk) = 0;
               end
            end
            if ismember(xF',O','rows')
               nS = nS + 1;
            end
            xD = x_LL + x.*x_DD;
            sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - sn_f(xD,snProb);
         end
      end
      if PLOTFLAG
         figure('Name','Difference RBF Surface vs. True Function');
         mesh(xGrid,xGrid,reshape(sFx,41,41))
      end

   elseif d == 3

      xGrid = linspace(0,1,21);
      N = length(xGrid)^d;
      sFx = zeros(N,1);
      fI  = true(N,1);        %Feasibility Index
      kk = 0;
      for ii = xGrid
         for jj = xGrid
            for ij = xGrid
               kk = kk+1;
               x = [ii jj ij]';
               xF = x_L + x.*x_D;
               if dLin + dCon(1) > 0
                  cVal = 0;
                  if dLin > 0
                     L = Prob.A*xF;
                     cVal = cVal+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
                  end
                  if dCon(1) > 0
                     C = cgo_c(xF, Prob, varargin{:});
                     nCon = nCon + 1;
                     cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
                  end
                  if cVal > 0
                     fI(kk) = 0;
                  end
               end
               if ismember(xF',O','rows')
                  nS = nS + 1;
               end
               xD = x_LL + x.*x_DD;
               sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - sn_f(xD,snProb);
            end
         end
      end

   elseif d == 4

      xGrid = linspace(0,1,11);
      N = length(xGrid)^d;
      sFx = zeros(N,1);
      fI  = true(N,1);        %Feasibility Index
      kk = 0;
      for ii = xGrid
         for jj = xGrid
            for ij = xGrid
               for ji = xGrid
                  kk = kk+1;
                  x = [ii jj ij ji]';
                  xF = x_L + x.*x_D;
                  if dLin + dCon(1) > 0
                     cVal = 0;
                     if dLin > 0
                        L = Prob.A*xF;
                        cVal = cVal+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
                     end
                     if dCon(1) > 0
                        C = cgo_c(xF, Prob, varargin{:});
                        nCon = nCon + 1;
                        cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
                     end
                     if cVal > 0
                        fI(kk) = 0;
                     end
                  end
                  if ismember(xF',O','rows')
                     nS = nS + 1;
                  end
                  xD = x_LL + x.*x_DD;
                  sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - sn_f(xD,snProb);
               end
            end
         end
      end

   end
   RMS = norm(sFx)/N;
   MASurfErr = max(abs(sFx));
   if all(fI)
      nF          = [];
      RMS_F       = RMS;
      MASurfErr_F = MASurfErr;
   else
      nF          = sum(fI);
      RMS_F       = norm(sFx(fI))/nF;
      MASurfErr_F = max(abs(sFx(fI)));
   end
else
   N           = [];
   nS          = [];
   nF          = [];
   RMS         = [];
   MASurfErr   = [];
   RMS_F       = [];
   MASurfErr_F = [];
end



% SAVE RESULTS

fMinIdx  = fIdx(1);
rngState = rand('state');
% Create struct output with warmstart information; save results on cgoSave.mat
Result.CGO.WarmStartInfo = saveCGO(2,'arbfMIP', Name,O,F,X,F_m,F00,Cc,...
                          nInit,Fpen,fMinIdx, rngState);

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

if dCon(1) > 0 & dCon(2) == 0
   cMin            = cgo_c(Result.x_k, Prob, varargin{:});
   nCon            = nCon + size(cMin,2);
   Result.c_k      = cMin;     % Constraint value at best x_k
   Result.ConstrEv = nCon;
   %for i=1:size(Result.x_k,2)
   %   cMin = [cMin,nlp_c(Result.x_k(:,i), Prob, varargin{:})];
   %   nCon = nCon + 1;
   %end
elseif dCon(1) > 0 & dCon(2) > 0
   cMin            = cgo_c(Result.x_k, Prob, varargin{:});
   nCon            = nCon + size(cMin,2);
   Result.CGO.c_k  = cMin;     % Noncostly constraint value at best x_k
   Result.CGO.nCon = nCon;
   Result.c_k      = CcMin;     % Costly constraint value at best x_k
   Result.ConstrEv = nFunc;
   %for i=1:size(Result.x_k,2)
   %   cMin = [cMin,nlp_c(Result.x_k(:,i), snProb, varargin{:})];
   %   nCon = nCon + 1;
   %end
elseif dCon(2) > 0
   Result.c_k        = CcMin;     % Costly constraint value at best x_k
   Result.ConstrEv   = nFunc;
else
   Result.ConstrEv   = 0;
   Result.c_k        = [];
end
%cMin = [];
%if dCon > 0
%   for i=1:size(Result.x_k,2)
%      cMin = [cMin,nlp_c(Result.x_k(:,i), Prob, varargin{:})];
%      nCon = nCon + 1;
%   end
%end
%if dLin > 0 & size(Result.x_k,2)==1
%   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at best x_k
%end
if dLin > 0
   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at all x_k
end
Result.CGO.snProb = snProb;
Result.CGO.gnProb = gnProb;
Result.CGO.fGoal  = fGoal;
Result.CGO.Its    = Its;
Result.Iter       = Iter;     % Number of iterations
Result.FuncEv     = nFunc;
Result.GradEv     = -1;
Result.HessEv     = -1;
Result.ConJacEv   = -1;
Result.ExitFlag   = 0;
Result.ExitText   = ['Tried ' num2str(nFunc) ' f(x), using ' ...
   num2str(n) ', startup ' num2str(nInit), ' ',ExDText];
Result.SolverAlgorithm = [Result.SolverAlgorithm ...
   '. Global solver ' globalSolver  ...
   '. Local solver ' localSolver];
if cpumax & convflag == 0
   Result.Inform   = 9;
   Result.ExitText = [Result.ExitText '. Max CPU reached. '];
   %elseif ~progress & convflag == 0
   %   Result.Inform   = 8;
   %   Result.ExitText = [Result.ExitText '. No progress for ' ...
   %                      num2str(nFunc-FLOWER) ' function evaluations'];
elseif convflag == 7
   Result.Inform   = convflag;
   Result.ExitText = [Result.ExitText, '. All feasible integers tried'];
else
   Result.Inform   = convflag;
end
Result.DIGIT       = TESTDIG;

%NHQ Save ExD test parameters
Result.ExD.xDIGITS     = xDIGITS;
Result.ExD.minXdist    = minXdist;
Result.ExD.RMS         = RMS;
Result.ExD.MASurfErr   = MASurfErr;
Result.ExD.N           = N;
Result.ExD.nS          = nS;
Result.ExD.RMS_F       = RMS_F;
Result.ExD.MASurfErr_F = MASurfErr_F;
Result.ExD.nF          = nF;

Result             = endSolve(Prob,Result);

if(usecgolib == 1)
  cgolib(101, Prob.CGOLIB.obs);
else
  tomsol(25) % Deallocates memory
end

% -------------- End of main ARBFMIP routine -------------------
% --------------------------------------------------------------
% --------------------------------------------------------------

% --------------------------------------------------------------
function snProbD1 = DefsnProbD1(snProb)
% --------------------------------------------------------------
% Create structure for optimization with added distance constraints
snProbD1                     = snProb;
snProbD1.ConsDiff            = 0;
dCon                         = snProb.mNonLin;
SCALE                        = snProb.SCALE;
snProbD1.mNonLin             = dCon+1;
snProbD1.c_L(dCon+1)         = 0;
snProbD1.c_U(dCon+1)         = inf;

if ~(dCon > 0 & SCALE > 0)
   snProbD1.FUNCS.c      = 'rbf_c';
   snProbD1.FUNCS.dc     = 'rbf_dc';
   snProbD1.FUNCS.d2c    = 'rbf_d2c';
end
if dCon == 0
   snProbD1.cNargin     = 0;
   snProbD1.dcNargin    = 0;
   snProbD1.d2cNargin   = 0;
elseif SCALE == 0
   snProbD1.c           = snProb.FUNCS.c;
   snProbD1.dc          = snProb.FUNCS.dc;
   snProbD1.d2c         = snProb.FUNCS.d2c;
   snProbD1.cNargin     = xnargin(snProbD1.c);
   if isempty(snProbD1.dc)
      snProbD1.dcNargin    = 0;
   else
      snProbD1.dcNargin    = xnargin(snProbD1.dc);
   end
   if isempty(snProbD1.d2c)
      snProbD1.d2cNargin   = 0;
   else
      snProbD1.d2cNargin   = xnargin(snProbD1.d2c);
   end
end
%----------------------------------------------------------------------------
function clusters=AdCluster(C,Prob)
%----------------------------------------------------------------------------
% Adaptive Clustering algorithm
%   maxLocalTry Maximal number of local searches from cluster points
%               If <= 0, glcCluster stops after clustering. Default 30
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

PriLev = 3;

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

if isfield(Prob.GO,'maxDistMin')
   maxDistMin   = Prob.GO.maxDistMin;
else
   maxDistMin = [];
end
if isempty(maxDistMin), maxDistMin = 0.01; end

% makeAllNorms returns the distance from any point
% to that points closest neighbour.
% minl1(i) = min (|| C(:,i) - C(:,j)||, j=1:m, j~=i, m=size(C,2))

minl1 = tomsol(29,C);     % minl1=makeAllNorms(C);
xprint(minl1,'minl1')

% Average distance any point has to its closest neighbour
meanDist=mean(minl1);
xprint(meanDist,'meanDist')
% meanDist will be an upper bound on maxDist

vIdx = minl1 >= meanDist;

% ixX is the set of points being considered for clustering
ixX=find(vIdx == 0);

% Points that are to be clustered will no longer be considered for trisection
% Now, for this subset, find the closest point to all point

minl = tomsol(29,C(:,ixX));

% Perhaps this min number should be changed so that
%         if max(minl) < 0.01 maxDist = 0.05
%         if  0.01 < max(minl) < 0.1 maxDist = 0.1
%         else maxDist = max(minl)

maxDist = max([max(minl) maxDistMin]);
% Avoid square root: maxDist = max([max(minl) 0.05^2]);

clusters     = cluster(C(:,ixX),maxDist);
clusterPnts  = find(clusters==0);
nPnt         = 1+length(clusterPnts);

if PriLev > 2
   fprintf('                      ')
   fprintf('maxDist %8.5f ',maxDist)
   fprintf('nPnt %3d ',nPnt)
   fprintf('pnts Checked %d ',length(ixX))
   fprintf('Total pnts %d ',size(C,2))
   fprintf('\n')
end

if nPnt > maxLocalTry
   while nPnt > maxLocalTry & maxDist < meanDist
      % Increase maxDist
      maxDist     = incDist*maxDist;
      clusters    = cluster(C(:,ixX),maxDist);
      clusterPnts = find(clusters==0);
      nPnt        = 1+length(clusterPnts);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(ixX))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
elseif nPnt < minLocalTry
   maxDist = decDist*maxDist;
   while nPnt < minLocalTry & maxDist > maxDistMin
      clustersSave    = clusters;
      clusterPntsSave = clusterPnts;
      clusters        = cluster(C(:,ixX),maxDist);
      clusterPnts     = find(clusters==0);
      nPnt            = 1+length(clusterPnts);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(ixX))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
      maxDist         = decDist*maxDist;
   end
   if nPnt > maxLocalTry
      % Too many clusters, use previous clustering
      clusters    = clustersSave;
      clusterPnts = clusterPntsSave;
      maxDist     = maxDist/decDist;
      nPnt        = 1+length(clusterPnts);
      if PriLev > 2
         fprintf('                      ')
         fprintf('Use previous clustering, nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(ixX))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
   maxDist        = maxDist/decDist;
elseif nPnt < maxLocalTry & incLocalTry
   maxDist = decDist*maxDist;
   while nPnt < maxLocalTry & maxDist > maxDistMin
      clustersSave    = clusters;
      clusterPntsSave = clusterPnts;
      clusters        = cluster(C(:,ixX),maxDist);
      clusterPnts     = find(clusters==0);
      nPnt            = 1+length(clusterPnts);
      if PriLev > 2
         fprintf('                      ')
         fprintf('maxDist %8.5f ',maxDist)
         fprintf('nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(ixX))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
      maxDist         = decDist*maxDist;
   end
   if nPnt > maxLocalTry
      % Too many clusters, use previous clustering
      clusters    = clustersSave;
      clusterPnts = clusterPntsSave;
      maxDist     = maxDist/decDist;
      nPnt        = 1+length(clusterPnts);
      if PriLev > 2
         fprintf('                      ')
         fprintf('Use previous clustering, nPnt %3d ',nPnt)
         fprintf('pnts Checked %d ',length(ixX))
         fprintf('Total pnts %d ',size(C,2))
         fprintf('\n')
      end
   end
   maxDist        = maxDist/decDist;
end


%----------------------------------------------------------------------------
function clusters=clusterX(X,maxDist)
%----------------------------------------------------------------------------

clusters=[];

[md idx]=min(sqrt(sum(X.^2,1)));
%[md idx]=min(sum(X.^2,1));

idx=idx(1); % We are only interested in one of these points.
% If the other points are in the same cluster they will be found any way.
% If they should start there own clusters, they will be found later on.

clusters(1) = idx;
clCnt       = 1;
ok          = 1;
k           = idx;

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
function GrpPrint(s,j,xC,Cand,GonBmin,GonBmax,Gsize,fx,doX)
%----------------------------------------------------------------------------
if ~ischar(s), s = num2str(s); end
fprintf('%s: Group %d ',s,j);
fprintf('onBmin %2d ',GonBmin(j));
fprintf('onBmax %2d ',GonBmax(j));
fprintf('Sz %2d ',Gsize(j));
fprintf('%6s  %2d ',Cand,xC);
if ~isempty(fx)
   fprintf('fx %14.10f ',fx(xC));
end
if size(doX,2) > 1
   fprintf('doX %9.5f/%d ',doX(xC,:));
else
   fprintf('doX %9.5f ',doX(xC));
end
fprintf('\n');

%----------------------------------------------------------------------------
function LicOK = blahaha(Solver)
LicOK = 1;

%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
function fnStarOpt(s,j,xC,Cand,GonBmin,GonBmax,Gsize,fx,doX)
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
function D = distx1x2(x1,x2,IntVars,Reals)
%----------------------------------------------------------------------------
if isempty(IntVars)
   D     = tomsol(30,x1,x2);
   D(2)  = 0;
else
   D(1)  = tomsol(30,x1(Reals),x2(Reals));
   D(2)  = sum(x1(IntVars)~=x2(IntVars));
end

%----------------------------------------------------------------------------
function fV = fnSStat(fnStar,xTry,fTry,gnProb,x_min,min_sn_y,X,epsX,IntVars,Reals)
%----------------------------------------------------------------------------
% Columns in fV, quantities computed for xTry
% 1: fnStar 2: onB(xTry) 3: s_n(xTry) 4: gn_f(xTry) 5: my(xTry)
% 6: doX 7: doM 8: doS if ~isempty(IntVars): 9: doX(2) 10: doM(2) 11: doS(2)

% Distance between new point and best point found or closest point in X

% doS  Distance from min_sn_y to new point xTry
doS             = distx1x2(min_sn_y,xTry,IntVars,Reals);

[onB, doX, doM] = statGN(gnProb.x_L,gnProb.x_U,xTry,x_min,X,epsX,IntVars,Reals);

if(gnProb.CGOLIB.usecgolib == 1)
  snNew = cgolib(303, gnProb.CGOLIB.rbf, xTry);
else
  snNew = tomsol(21, xTry);
end

if isinf(fnStar)
   my = -1/fTry;
else
   z  = (snNew-fnStar)^2;
   if z~=0 & fTry ~= 0
      my = -1/(fTry*z);
   else
      my = inf;
   end
end
if isempty(IntVars)
   fV = [fnStar,onB,snNew,fTry,my,doX(1),doM(1),doS(1)];
else
   fV = [fnStar,onB,snNew,fTry,my,doX(1),doM(1),doS(1),...
      doX(2),doM(2),doS(2)];
end

%----------------------------------------------------------------------------
function fx = fTrue(xTry,x_L,x_D,Prob, varargin)
%----------------------------------------------------------------------------
% New point in original space
if Prob.CGO.SCALE
   O_new   = tomsol(9, x_L, xTry, x_D);
else
   O_new   = xTry;
end
% HKH - No use of costly Cc(x) in CcNew right now
[fx,CcNew] = cgo_fc(O_new, Prob, varargin{:});

%% global NLP_x NLP_f NARG
%NLP_x=[]; NLP_f=[]; NARG = [];
%fx  = nlp_f(O_new, Prob);

%----------------------------------------------------------------------------
function fnSPrint(XJ,fJ,min_sn,DebugARBF)
%----------------------------------------------------------------------------
nTry = size(fJ,1);
for i = 1:nTry
   fprintf('%2d: ',i);
   fnStar = fJ(i,1);
   fprintf('fnS %10.6f ',fnStar);
   onB = fJ(i,2);
   fprintf('[%d] ',onB);
   if DebugARBF
      UP = 0;
      if i > 1
         if fJ(i,end) > fJ(i-1,end)
            UP = 1;
         end
      end
      if UP == 0
         fprintf('fx %13.8f ',fJ(i,end));
      else
         fprintf('fX %13.8f ',fJ(i,end));
      end
   end
   snNew = fJ(i,3);
   fprintf('sn %11.6g ',snNew);
   % fprintf('Le %7.3f ',L1);
   % fprintf('L %7.3f ',L2);
   if size(fJ,2) < 12
      fprintf('doM %6.3f ',fJ(i,7));
      fprintf('doX %6.3f ',fJ(i,6));
      fprintf('doS %6.3f ',fJ(i,8));
   else
      fprintf('doM %6.3f/%d ',fJ(i,7),fJ(i,10))
      fprintf('doX %6.3f/%d ',fJ(i,6),fJ(i,9))
      fprintf('doS %6.3f/%d ',fJ(i,8),fJ(i,11))
   end
   my = fJ(i,5);
   fprintf('my %10.4f ',my);
   if my > 0
      myLog = log10(my);
   elseif my == 0
      myLog = inf;
   else
      myLog = log10(abs(my));
   end
   fprintf('myL %5.2f ',myLog);
   fTry = fJ(i,4);
   gnLog = log10(max(1E-300,min_sn-fTry));
   if gnLog <= -300
      fprintf('gnL %6.2f ',0);
   else
      fprintf('gnL %6.2f ',gnLog);
   end
   fprintf('gn_f %11.6f ',fTry);
   xTry = XJ(:,i);
   xprint(xTry);
end

%----------------------------------------------------------------------------
function D = normDist(x1,x2,dSQRT,IntVars,Reals)
%----------------------------------------------------------------------------
if isempty(IntVars)
   D(1)  = tomsol(30,x1,x2)/dSQRT;
   D(2)  = 0;
else
   D(1)  = tomsol(30,x1(Reals),x2(Reals))/dSQRT;
   D(2)  = sum(x1(IntVars)~=x2(IntVars));
end

% --------------------------------------------------------------
function [idea, rbfType, N, infStep, fStarRule, ...
   DeltaRule, MaxCycle, eps_sn, TargetMin, AddSurfMin, ...
   usecgolib] = getARBFProbVars(Prob, x_L, x_U)
% --------------------------------------------------------------


if isempty(Prob.CGO)
   %idea         = [];
   rbfType      = []; N         = []; infStep    = [];
   fStarRule    = []; DeltaRule    = []; MaxCycle  = []; eps_sn     = [];
   TargetMin    = []; AddSurfMin   = [];
else
   %if isfield(Prob.CGO,'idea')
   %   idea = Prob.CGO.idea;
   %else
   %   idea = [];
   %end
   if isfield(Prob.CGO,'rbfType')
      rbfType = Prob.CGO.rbfType;
   else
      rbfType = [];
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
   if isfield(Prob.CGO,'fStarRule')
      fStarRule = Prob.CGO.fStarRule;
   else
      fStarRule = [];
   end
   if isfield(Prob.CGO,'DeltaRule')
      DeltaRule = Prob.CGO.DeltaRule;
   else
      DeltaRule = [];
   end
   if isfield(Prob.CGO,'MaxCycle')
      MaxCycle = Prob.CGO.MaxCycle;
   else
      MaxCycle = [];
   end
   if isfield(Prob.CGO,'eps_sn')
      eps_sn = Prob.CGO.eps_sn;
   else
      eps_sn = [];
   end
   if isfield(Prob.CGO,'TargetMin')
      TargetMin = Prob.CGO.TargetMin;
   else
      TargetMin = [];
   end
   if isfield(Prob.CGO,'AddSurfMin')
      AddSurfMin = Prob.CGO.AddSurfMin;
   else
      AddSurfMin = [];
   end
   if isfield(Prob.CGO, 'usecgolib')
     usecgolib = Prob.CGO.usecgolib;
   else
     usecgolib = 1;
   end
end

% if isempty(idea),         idea = 1; end
% HKH always set idea = 1
idea = 1;
if isempty(rbfType),      rbfType = 2; end
if isempty(infStep),      infStep = 0; end
if isempty(fStarRule),    fStarRule = 1; end
if isempty(DeltaRule),    DeltaRule = 1; end
if isempty(MaxCycle),     MaxCycle = 10; end
if isempty(eps_sn),       eps_sn    = 1E-7; end
if isempty(AddSurfMin),   AddSurfMin = 0; end
if isempty(TargetMin),    TargetMin = 3; end

if fStarRule == 3
   infStep = 1; % Must have infStep in strategy 3, (Gutmann I)
   if isempty(N), N = 1; end
else
   if isempty(N), N = 2; end
end


%
% 061201  hkh  New arbfMIP algorithm implemented in 1st version
% 070308  hkh  Major revision to handle MINLP
% 070308  hkh  Define backupSolver1 and backupSolver2, safe snSolve computation
% 071006  hkh  Do init of random generator only if no warm start, NaN no init
% 071007  hkh  Save random generator state in rngState, reinit rand if WarmStart
% 071007  hkh  Add TargetMin parameter, 4 strategies to select minimum
% 071010  hkh  Major revision of ARBF algorithm
% 080110  hkh  Revision of comments and output, streamline with ego and rbfSolve
% 080410  hkh  Complete revision in modules, new experimental design
% 080414  hkh  Prob.GO input in Result.CGO, change CGOGLobalProb call
% 080414  hkh  Add nTrial, CLHMethod to Result.CGO
% 080417  hkh  Use relative constraint violation in fPen computation
% 080420  hkh  Change fnStar strategy, use fnS=min_sn-vJ*|min_sn|, with larger grid vJ
% 080420  hkh  Use 1 or more additional fnStar values when adding extra values in the grid
% 080421  hkh  If FASTCODE=1, avoid optimize for unnecessary fnStar values down to -inf
% 080507  hkh  Add modN == -11, adding infStep before min_sn in global grid
% 080617  frhe Made REPLACE>1 transformation monotonic
% 080617  frhe log10 used in REPLACE>1 in accordance with help
% 080619  frhe CGOLIB calls added
% 080630  hkh  Use cgo_fc instead of nlp_f for costly f(x) (and costly c(x))
% 080711  hkh  Use routine saveCGO for warm start info saving for all solvers
% 080711  hkh  Save warm start info every iteration in cgoSave.mat (cgoSave1)
% 080711  hkh  Avoid estimation of Hessian in endSolve setting HessEv = -1
% 081104  hkh  Avoid crash in eig(H) when sn_f/sn_g/sn_H is bad with NaN values
% 081105  hkh  Make getProbVars solver dependent
% 081105  hkh  Make all extra output dependent on PriLev > 1
% 081105  hkh  Add input/output nTrial in call to getProbVars
% 090813  med  repmat call corrected
% 090824  hkh  Minor mlint revision
% 090918  hkh  Send Fpen in snProb.CGO.F to snSolve
% 091020  hkh  Use x_k(:,1) from special call to multiMin

%HKH::::
% Ev. Stoppa om alla punkter utom 1: är == inf-punkten? likely minimum

% Kolla Red=(predicted-actual)/actual för alla punkter. Om globalt steg ger alla punkter 
% har Red < tol (1E-4, eps_f?) så stoppa. Speciellt inf-min skall ha denna noggrannhet

% Måste hålla koll på om lokalt min nåtts, och då starta speciell global-fas
% Använda ExpI för RBF? inf-min? flera lokala min från inf-min direkt?
% Skall man köra linjära bivillkor för att hålla undan nuvarande lokalt min man funnit,
% eller kanske addera en barriär eller penalty för denna

% Hålla en lista av de minimum man hittat

% Kanske kolla på RBF-ytan, lägga ett bivillkor runt min hittat så att nytt min_sn_y nås
% någonstans på annat ställe

% Implementera ALLA RBF-funktioner!

% ev. log av data med stora skalor. jmf Goldstein med o utan log
% Transformationer a la ego?

% Visa testresultat a la shoemakers doktorand

% Går det att använda storlek på my_n, dvs min för inf-problemet?

