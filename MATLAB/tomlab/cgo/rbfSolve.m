% rbfSolve is based on the RBF algorithms presented in
% 1. Bjorkman, Holmstrom: Global Optimization of Costly Nonconvext Functions
%    Using Radial Basis Functions, Optimization and Engineering 1,373-397, 2000
% 2. Hans-Martin Gutmann: A radial basis function method for global
%    optimization, Journal of Global Optimization, 19,201:207, 2001.
% 3. Hans-Martin Gutmann: Radial Basis Function Methods for Global Optimization,
%    Ph.D. Thesis, Cambridge University, Cambridge, UK, September 2001.
%
% rbfSolve solves problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Some or all x may be integer valued as specified by other input variables
%
% f(x) are assumed to be a costly function
% c(x) are assumed to be cheaply computed
%
% Any set {j} of costly constraints can be treated by adding penalty terms to 
% the objective function in the following way:
%      p(x) = f(x) + SUM_j w_j * max (0,c^j(x)-c_U^j, c_L^j-c^j(x)),
% where weighting parameters w_j > 0 have been added.
% The user then returns p(x) instead of f(x) to the CGO solver.
%
% Calling syntax:
%
% function Result = rbfSolve(Prob, varargin)
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
%   WarmStart If true, >0, rbfSolve reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%             If Prob.CGO.WarmStartInfo has been defined through a call to
%             WarmDefGLOBAL, this field is used instead of the cgoSave.mat file.
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   user      User field used to send information to low-level functions
%   PriLevOpt Print Level
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%   PriLevSub Print Level in subproblem solvers, see help in snSolve, gnSolve
%   f_Low     Lower bound on the optimal function value. If defined, used to
%             restrict the target values into interval [f_Low,min(surface)]
% -----------------------
% Fields in Prob.optParam
% -----------------------
%             Defines optimization parameters. Fields used:
%  MaxFunc    Maximal number of costly function evaluations. Default 300
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then set MaxFunc = MaxFunc + nFunc
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
%             except for pure IP problems, then max(GO.MaxFunc, MaxIter);
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
% AddMP        If = 1, add the midpoint as extra point in the corner strategies
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
% SCALE        0-Original search space (Default if any integer values)
%              1-Transform search space to unit cube (Default if no integers).
%
% REPLACE      0-No replacement (Default for constrained problems)
%              1-Large function values are replaced by the median
%              >1 - Large values Z are replaced by new values
%              Replacement: Z:= FMAX + log10(Z-FMAX+1), where
%              FMAX = 10^REPLACE, if min(F) < 0
%              FMAX = 10^(ceil(log10(min(F)))+REPLACE), if min(F) >= 0
%              (Default REPLACE = 5 if no linear or nonlinear constraints)
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
%              --- Special RBF algorithm parameters
%
% rbfType      Type of radial basis function.
%              1 - Thin Plate Spline, 2 - Cubic (Default).
%              3 - Multiquadric,      4 - Inverse multiquadric
%              5 - Gaussian,          6 - Linear
%
% idea         Type of search strategy on the surface, idea 1 (cycle of N+1=5
%              points as default, with different target values fnStar) or
%              idea 2 (cycle of 4 points in alpha, which implicitly gives
%              the target value fStar).
%              Default is idea 1, cycle length N = Prob.CGO.N.
%
% N            Cycle length in idea 1 (Default N=1 for fStarRule 3, otherwise N=4)
%              or in idea 2 (always N=3)
%
% infStep      If =1, add search step with target value -inf first in cycle
%              Default 0. Always =1 for the case idea=1, fStarRule=3
%
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
% DeltaRule    1 = Skip large f(x) when computing f(x) interval Delta
%              0 = Use all points. Default 1.
%
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
%
% eps_sn       Relative tolerance used to test if minimum of surface, min_sn, is
%              sufficiently lower than the best point found (fMin).  Default = 1E-7
% MaxCycle     Max number of cycles without progress before stopping. Default 10
% NOISE        0 = Objective function is not noisy. 
%              1 = Objective is noisy, use an approximation of
%                  of the surrogate model.
% usecgolib    1 = Use CGOLIB
%              0 = Do not use CGOLIB
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
% ------------------
% Fields in Prob.MIP
% ------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:
%             If empty, all variables are assumed non-integer
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% --------------------------------------------------------------------
% varargin    Additional parameters to rbfSolve are sent to the costly f(x)
% --------------------------------------------------------------------
%
% -----------------
% OUTPUT PARAMETERS
% -----------------
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
%           4 = No new point sampled for N iteration steps
%           5 = All sample points same as the best point for N last iterations
%           6 = All sample points same as previous point for N last iterations
%           7 = All feasible integers tried
%           8 = No progress for MaxCycle*(N+1)+1 function evaluations
%               (>MaxCycle cycles, input CGO.MaxCycle)
%           9 = Max CPU Time reached
%
% To make a warm start possible, rbfSolve saves the following information in
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
% Cycle     Cycle steps global to local. infStep is marked -1,
%           0 to N-1 are global steps. Last step N in cycle is surface minimum
% R         If the letter R is printed, the current step is a RESCUE step, i.e.
%           the new point is already sampled in a previous step, instead the
%           surface minimum is used as a rescue
% fnStar    Target value fn_star (alpha for idea 2).
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
% xOpt      User-given global optimum Prob.x_opt (if defined
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
%
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
% --- Row 4 (if PriLev > 3)
% -----------
% snNew-min_sn
% snNew-fnStar
% snNew-fNew
% myNew
% fRed
% hn
% hnErr
% LipU      Maximum Lipschitz constant for initial set X
% LipL      Minimum Lipschitz constant for initial set X
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
%    Prob.optParam.MaxFunc = 400;  % Change max number of function evaluations
%    Prob.GO.MaxFunc       = 20000;% 20000 global function values in subproblem
%
% Direct solver call:
%    Result = rbfSolve(Prob);
%    PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('rbfSolve',Prob,2);
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
% To make a restart, just set the restart flag, and call rbfSolve once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('rbfSolve',Prob,2);
%
% Another option for warm start is to put the restart information from the
% Result structure to the Prob structure and then call a CGO solver again:
%    Prob = WarmDefGLOBAL('rbfSolve',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('rbfSolve',Prob,2);
%
% It is also possible to run another CGO solver in the warm start:
%    Prob = WarmDefGLOBAL('rbfSolve',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('ego',Prob,2);
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 12, 2000. Last modified Oct 20, 2009.

function Result = rbfSolve(Prob, varargin)

if nargin < 1
   error('rbfSolve needs one parameter, the structure Prob');
end
global NLP_x NLP_f NARG

DEBUG       = 0;
DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt
SolveInf    = 0;  % If to solve the inf problem first in every step, time consuming

time                   = fix(clock);
solvType               = checkType('glc');
Prob                   = ProbCheck(Prob,'rbfSolve',solvType);
Prob                   = iniSolve(Prob,solvType,0,0);

Result                 = ResultDef(Prob);
Result.Solver          = 'rbfSolve';
Result.SolverAlgorithm = 'Radial Basis Function Interpolation';

% Pick up bounds from Prob structure and check if OK
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds

if isempty(x_L) | isempty(x_U)
   disp('rbfSolve requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'rbfSolve requires both lower and upper variable bounds';
   Result.DIGIT    = [];
   Result.CGO.Its  = [];
   Result          = endSolve(Prob,Result);
   return;
end
if ~(all(isfinite(x_L)) & all(isfinite(Prob.x_U)))
   disp('rbfSolve only solves box bounded problems, where both');
   disp('lower bounds Prob.x_L and upper bounds Prob.x_U are finite');
   Result.ExitFlag = 1;
   Result.ExitText = 'Problem not box bounded, variable bounds on x not finite';
   Result.DIGIT    = [];
   Result.CGO.Its  = [];
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
PEN=Prob.CGO.PEN;
SAVE = 1;
% ------------------------------------ ----------------------
% Set possibly redefined GOlocalSolver and DIRECT into Prob.GO
% ----------------------------------------------------------
Prob.GO.localSolver  = GOlocalSolver;
Prob.GO.DIRECT       = DIRECT;

% -------------------------------------------------
% Get special RBF variables from Prob.CGO structure
% -------------------------------------------------
[idea, rbfType, N, infStep, fStarRule, ...
   DeltaRule, MaxCycle, eps_sn, TargetMin, AddSurfMin, ...
   usecgolib, NOISE] = getRBFProbVars(Prob, IntVars);

% -------------------------------------------------------------
% Get variables from Prob structure
% Redefine GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3
% and set back in Prob.GO
% -------------------------------------------------------------

[Prob ,MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, ... 
          fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX, x_LL, x_UU, ...
          x_DD, pDist, dLin, dCon, Percent, nSample, nTrial, AddMP, REPLACE,...
          backupSolver1, backupSolver2, xOptS ] = ...
          getProbVars(Prob, x_L,x_U,x_D, SCALE,globalSolver,IntVars, ...
          Percent, nSample, nTrial, AddMP, REPLACE, GOMaxFunc, GOMaxIter, ...
          maxFunc1, maxFunc2, maxFunc3,'rbfSolve');

%  STEP 1, Initialization


nmax      = [];       % Only used in idea 1, in rbfSolve

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
   [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, nCon, nSample,...
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

Prob.CGOLIB = struct('usecgolib', usecgolib);

% CGOLIB Guide:
%
% Man skapar först ett observation set. Ett observation set består av
% punkter (kallas rader) och funktionsvärden (kallas kolumner). Man kan
% alltså ha flera funktionsvärden för samma punkt, och det används till
% REPLACE-transformer och log-transformer bl.a.
%
% Ex.: Vi har punkter X (lagrade kolonnvis) och funktionsvärden F (en stor
% kolonnvektor med lika månge element som X har kolonner, d.v.s. ett
% funktionsvärde per punkt):
%
% setid = cgolib(100, size(X,1), size(X,2), 1);
% Detta skapar ett observation set som till en början får plats med
% size(X,2) rader (punkter) med dimension size(X,1) och 1 kolonn (funktionsvärden).
% Detta betyder inte att några rader eller kolonner har skapats, bara att
% det har allokerats plats för att lagra så många. Detta utvidgas
% automatiskt om man lägger till fler rader eller kolonner.
%
% setid är en identifier till denna mängd.
%
%
% colnr = cgolib(103, setid, 1)
% Detta _skapar_ en datakolonn. (Man skiljer på datakolonner och transformationskolonner).
% 
% colnr är kolumnens nummer.
%
%
% cgolib(102, setid, X, colnr, F);
% Detta fyller mängden med punkter (X) och funktionsvärden (F). F skulle
% kunna vara värdena från bivillkoren. Man skulle kunna ha flera kolonner,
% en för varje bivillkor.
%
%
% cgolib(199, setid)
% Detta skriver ut hela observation settet på skärmen. Är inte bra i
% färdiga programmet, men kan användas för att debugga.
%
%
% Sedan kan vi skapa ytor och koppla till en kolumn i ett obs. set:
% rbfid = cgolib(300, setid, colnr)
% Dett skapar en RBF-yta från kolumn colnr i set setid. rbfif är identifier
% till den nya ytan.
%
%
% Man kan sätta options med:
% rbfoptions.PHIFUNCTION=3; % välj multiquadric
% cgolib(350, rbfid, rbfoptions);
%
%
% Sedan måste man interpolera ytan innan den används:
% cgolib(302, rbfid);
%
%
% Man kan beräkna funktioner: sn_f, sn_g, gn_g, gn_f, h1n_f, m.fl. genom
% att göra olika anrop till cgolib. Se hjälpen (cgolib_mex.c) för alla
% alternativ och argument.
%
%
% Man kan checka interpolationen genom att köra:
% cgolib(323, rbfid); 
% Returnerar 0 om ok, -1 om inte ok, -2 om inte interpolerad yta.
%
%
% Om man lägger till en punkt i observation set så kommer ytan att meddelas
% att settet som den bygger på ändrats, och man måste interpolera om den
% för att man ska kunna räkna ut sn_f, m.fl. igen.
%
%
% Man kan skapa transformationskolumner också, som man sedan kan basera
% RBF-ytor på. Transformationskolumner är transformer av andra kolumner. En
% transformationskolumn kan vara en transformation av en annan
% transformationskolumn. Man måste uttryckligen räkna om
% transformationskolumnerna efter en ändring i settet. Det finns
% funktioner: cgolib(111, cgolib(112, som har med detta att göra.
%
%
% Det finns exempel på allt jag nämnt i koden nedan på olika ställen. Sök
% bara på 'cgolib('.

% -----------------------------------------
% Initial RBF interpolation surface
% -----------------------------------------
% TOMSOL INIT, send F_m to Fortran, not F
clear('tomsol');

if(usecgolib == 1)

  % Create the observation set
  InitSize = max(2*size(X,2), min(300,MaxFunc));
  Prob.CGOLIB.obs = cgolib(100, size(X,1), InitSize, 4);
  
  % Set option: MINDISTUPDATE
  obsoptions.MINDISTUPDATE = 1;
  cgolib(150, Prob.CGOLIB.obs, obsoptions);
  
  % Create the data column
  Prob.CGOLIB.datacol = cgolib(103, Prob.CGOLIB.obs, 1);
  
  % Add row data (original data)
  cgolib(102, Prob.CGOLIB.obs, X, Prob.CGOLIB.datacol, F);

  % To handle noise, add another column
  % See updateNoiseColumn for some comments on NOISE.
  if(NOISE == 1)
      Prob.CGOLIB.noisecol = cgolib(103, Prob.CGOLIB.obs, 1);
      Prob.CGOLIB.noiserbf = cgolib(300, Prob.CGOLIB.obs, Prob.CGOLIB.noisecol);    
      % List of etas to try
      etas = [0.001 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
      updateNoiseColumn(etas, Prob);
      
  else
      Prob.CGOLIB.noisecol = Prob.CGOLIB.datacol;
  end
  
  % If we are using the REPLACE transform, then create
  % a new column with this transform and base the RBF
  % surface on that column.
  
  % Transformation column doesn't exist yet.
  if(REPLACE == 1)
    Prob.CGOLIB.repcol = cgolib(111, Prob.CGOLIB.obs, 1, Prob.CGOLIB.noisecol);
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.repcol;
    cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
  elseif(REPLACE == 0)
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.noisecol;    
  elseif(REPLACE > 1)
    Prob.CGOLIB.repcol = cgolib(111, Prob.CGOLIB.obs, 6, ...
				Prob.CGOLIB.noisecol, ...
				{abs(REPLACE)});
    Prob.CGOLIB.rbfcol = Prob.CGOLIB.repcol;
    
%    Prob.CGOLIB.repcol = cgolib(111, Prob.CGOLIB.obs, 2, ...
%                                Prob.CGOLIB.datacol);
%   Prob.CGOLIB.rbfcol = Prob.CGOLIB.repcol;

    cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
else
    error('Invalid replace');
  end
  
  % Now create RBF surfaces
  % 1. Thin-plate
  % 2. Cubic
  % 3. Multiquadric
  % 4. Inverse multiquadric
  % 5. Gaussian
  % 6. Linear

  % Compute minimum distance within the set X
  mindist = cgolib(114, Prob.CGOLIB.obs);
  
  % HKH Could compute all 6 types if to use switching of RBF type
  % for i = 1:6
  for i = rbfType
    % Create a surface
    Prob.CGOLIB.rbftype(i) = cgolib(300, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
    
    % Choose RBF type
    rbfoptions = struct('');
    rbfoptions(i).PHIFUNCTION=i;
    rbfoptions(i).LUMETHOD=1;
    
    % Set some minimum distance dependent gamma to gain
    % numerical stability when solving for the interpolants later.
    if(i == 3 || i == 4)
        rbfoptions(i).PHIGAMMA = sqrt(mindist);
    elseif(i == 5)
        rbfoptions(i).PHIGAMMA = 1/mindist;
    end
    
    % Apply the options
    cgolib(350, Prob.CGOLIB.rbftype(i), rbfoptions(i));
  
    % Interpolate the new surface:
    cgolib(302, Prob.CGOLIB.rbftype(i));
    
    % Check the interpolations! Maybe we shouldn't return just because ONE
    % interplation fails. Maybe the failing interpolations should just be
    % removed.
    if(cgolib(323, Prob.CGOLIB.rbftype(i)) ~= 0)
        if(i == rbfType)
            cgolib(101, Prob.CGOLIB.obs); % Destroy the observation set and all children surfaces
            Result.ExitFlag = 1;
            Result.ExitText = 'Initial interpolation failed';
            Result.DIGIT    = [];
            Result.CGO.Its  = [];
	    fprintf('Rank of X matrix %d, should be %d\n',rank(X),d);

            % Create struct output with warmstart information; 
            % Save results on cgoSave.mat
            fMinIdx  = fIdx(1);
            rngState = rand('state');
            Result.CGO.WarmStartInfo = saveCGO(2,'rbfSolve',Name,O,F,X,F_m,...
                       F00,Cc,nInit,Fpen,fMinIdx, rngState);
            Result   = endSolve(Prob,Result);
            return
        else
            fprintf(['Warning: Interpolation of rbf surface type ' num2str(i) ' might be bad.\n']);
        end
    end
  end
  
  % For now, choose one of the RBF types as the one to use.
  % NOTE! Prob.CGOLIB.rbf should NOT contain the value according
  % to the list above, where 1 is thin-plate, 2 is cubic and so on.
  % Prob.CGOLIB.rbf is the RBF surface IDENTIFIER from CGOLIB. So it must
  % be taken from the Prob.CGOLIB.rbftype array.
  Prob.CGOLIB.rbf = Prob.CGOLIB.rbftype(rbfType);
  
  % The remaining part of rbfSolve is written for handling only one rbf
  % type, yet. This rbf type is given as the rbf surface identifier
  % through xProb.CGOLIB.rbf. The identifiers are now stored in
  % Prob.CGOLIB.rbftype(:). Note, that xProb is the problem that is solved
  % by a sub solver, for example snProb, or snProb.GO.ProbL. So before
  % solving a sub-problem, choose the rbfType by setting xProb.CGOLIB.rbf
  % and xProb.GO.ProbL.CGOLIB.rbf to the wanted rbf type identifier.
  
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
   snProb.optParam.IterPrint          = 1;
   snProb.PriLevOpt                   = DebugPriLev;
   snProb.GO.ProbL.optParam.IterPrint = 1;
   snProb.GO.ProbL.PriLevOpt          = DebugPriLev;
end
% --------------------------------------------------

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
   gnProb.optParam.IterPrint          = 1;
   gnProb.PriLevOpt                   = DebugPriLev;
   gnProb.GO.ProbL.optParam.IterPrint = 1;
   gnProb.GO.ProbL.PriLevOpt          = DebugPriLev;
end
gnProb.CGO.idea           = idea;
gnProb.GO.ProbL.CGO.idea  = idea;

%--------------------------------------------------------------

% HKH Not used yet
% dSQRT     = sqrt(d); % Used in distance calculation

% Set parameters used in gnSolve
gnProb.CGO.globalSolver = globalSolver;
gnProb.CGO.TargetMin    = TargetMin;
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

% Becase new sn_f,sn_g tests on CGO.fnStar
gnProb.CGO.fnStar = NaN;
if isempty(IntVars)
   snProb = ProbCheck(snProb,globalSolver,checkType('con'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('glc'));
else
   snProb = ProbCheck(snProb,localSolver,checkType('minlp'));
   gnProb = ProbCheck(gnProb,globalSolver,checkType('minlp'));
end
if 1
   % Improve speed
   snProb.FAST             = 1;
   snProb.GO.ProbL.FAST    = 1;
   gnProb.FAST             = 1;
   gnProb.GO.ProbL.FAST    = 1;
   % gnProb.GO.ProbL.NumDiff = 1;
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
% Find max Lipschitz constant in initial set

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
% nInit = Number of points used for startup, now or previous run
% Iter  = Number of iterations to obtain new sample points
% modN  = Cycle step

modN       = -1-infStep;
addlocStep = 0;

Iter       = 0;
Update     = 1;
fnStar     = Inf;
alpha     = NaN;
if idea == 2
   if dLin + dCon(1) > 0
      OldAlpha = 0; % If false, use a new Idea 2 alpha  strategy
   else
      OldAlpha = 1; % If false, use a new Idea 2 alpha  strategy
   end
end


FLOWER    = fIdx(1);
fMinIter  = 0;             % Only used in alpha computation
fDiff     = Inf;
fDiffOld  = NaN;
NOUPDATE  = 0;             % Used in test for Inform = 4
SAME1     = 0;             % Used in test for Inform = 5
SAME2     = 0;             % Used in test for Inform = 6
%HKH NOTMaxCycle  = 5;             % Used in tests for Inform 4,5,6
O_pre     = inf*ones(d,1);

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
   snOptH = nlp_H(xOptS,snProb);
   snOptE = sum(eig(snOptH)<0);
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
   fprintf('\n');
   maxF = max(F);
   fprintf('  max(F) %9.5g', maxF);
   fprintf(' med(F) %9.5g', median(F));
   fprintf(' rng(F) %9.5g', maxF-fMin);
   fprintf(' pDist %9.5g ', pDist);
   fprintf('LipU %f. ',LipUpp)
   fprintf('LipL %f. ',LipLow)
   fprintf('\n');
   xprint(O_min,'  xMin:',' %12.8f',8)
   if ~isempty(xOptS)
      if PriLev > 2
         xprint(Prob.x_opt,'  xOpt:',' %12.8f',8)
         fprintf('  SumXO (sum||xOpt-X||) %f ',SumXO);
         fprintf('MeanXO (mean||xOpt-X||) %f ',MeanXO);
         doO = min(tomsol(30,x_min,xOptS));
         fprintf(' doO (||xOpt-xMin|| %f ',doO);
         fprintf('\n');
         if PriLev > 3
            fprintf('  dXO (min||xOpt-X||) %f ',dXO);
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
   %fprintf('\n');
end

% -------------- Result saving -------------------

TIME0    = Prob.TIME0;
DIGIT    = 0.1;
TESTDIG  = [];
convflag = 0;
cpumax   = 0;
progress = 1;
FLOW     = nFunc;
%HKH NOT YET USED
if 0
   LocErr    = Inf;
   RelLocErr = Inf;
   PROGRESS  = Inf*ones(N+1,1);
   GlobErr(1:N)    = Inf;
   RelGlobErr(1:N) = Inf;
   GlobImp(1:N)    = Inf;
end

snProb.ixDist          = [];
snProb.GO.ProbL.ixDist = [];
gnProb.ixDist          = [];
gnProb.GO.ProbL.ixDist = [];


%NHQ 
% PLOTFLAG 0/1. Flag used to control plots. 
% If set to 1, plots are produced when d == 2.
PLOTFLAG = 0;
if d ~= 2
   PLOTFLAG = 0;
end
if PLOTFLAG
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

% Try loop and catch if stopped, copy cgoSave1/2 to cgoSave
while nFunc < MaxFunc & convflag == 0

   if nFunc-FLOW > MaxCycle*(N+1), progress = 0; break; end
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

   Iter                       = Iter + 1;

   % Set parameters used in snSolve
   snProb.snP                 = Iter;
   snProb.x_0                 = x_min;
   snProb.CGO                 = gnProb.CGO;
   snProb.CGO.globalSolver    = globalSolver;
   snProb.CGO.localSolver     = localSolver;
   snProb.CGO.fnStar          = NaN;
   % HKH modN not updated yet, make preliminary update for snSolve
   % Needed? not used in arbf
   snProb.CGO.modN            = mod(modN+1,N+1);
   snProb.PriLev              = PriSub;
   snProb.GO.ProbL.CGO        = snProb.CGO;
   
   snResult                   = snSolve(snProb);
   min_sn                     = snResult.f_k;


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
      disp('Min on sn-surface')
      disp(snResult.x_k)
      figure(rbfFIG)
      plot3(snResult.x_k(1),snResult.x_k(2),min_sn,'r*')

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


   if Feasible & min_sn-fMin > 1E-6*max(1,abs(fMin)) & NOISE == 0
      fprintf('WARNING min_sn %15.12f > fMin %15.12f. ',min_sn,fMin);
      if ~isempty(backupSolver1)
         fprintf('Rerun with %s\n',backupSolver1);
         snProb.CGO.globalSolver = backupSolver1;
         snProb.CGO.localSolver  = backupSolver1;
         if strcmpi(backupSolver1,'oqnlp')
            snProb.OQNLP.options.RANDOMSEED = 12345;
         end
         snProb.GO.ProbL.CGO     = snProb.CGO; % Just for safety, not needed
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
            snProb.GO.ProbL.CGO     = snProb.CGO; % Just for safety, not needed
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

   % min_snc and c0 NOT USED NOW
   %% Transform to original space
   %if SCALE
   %   min_snc = tomsol(9, x_L, min_sn_y, x_D);
   %   c0      = tomsol(9, x_L, x_min, x_D);
   %else
   %   min_snc = min_sn_y;
   %   c0      = x_min;
   %end


   % Choose a target value fnStar or compute alpha (implicitly gives fnStar)

   modN = modN + 1;
   if modN > N
      if FLOWER < nFunc - N
         % No improvement during the last cycle. Currently do nothing special
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
      [fnStar, ucOK, nmax] = getfnStar(modN, fnStar, n, nmax, nFunc, nInit, N,...
          DeltaRule, REPLACE, F, F_m, min_sn, fMin, ...
          f_Low, eps_sn, fStarRule);
      gnProb.CGO.fnStar = fnStar;
   elseif idea == 2
      [alpha, ucOK] = getAlpha(modN, alpha, F, F_m, min_sn, fMin, ...
          f_Low, eps_sn, Update, Iter, fMinIter, fDiffOld, OldAlpha);
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
         disp(xLoc)
      end
      % Remove the current interior point from the list of interior points in (xOL,fOL)
      xOL = xOL(:,2:end);
      fOL = fOL(2:end);
      % Not set (xInf,fInf), use same as last iteration
      %xInf = xLoc;
      %fInf = Inf;
   end
   if (modN == N) & ucOK
      dXsny               = min(tomsol(30,min_sn_y,X));
      if SolveInf
         NLP_x=[]; NLP_f=[]; NARG = [];
         gnProb.CGO.modN     = -1;
         gnProb.CGO.fnStar   = 0;
         gnProb.CGO.alpha    = 0;
         gnProb.PriLev       = PriSub;
         if dXsny < 1E-10 & length(IntVars) == d
            k                      = size(gnProb.CGO.X,2);
            XX                     = mipFeasible(gnProb, 4000, 1, 1, 0);
	    snResult.multiMin.xOpt = XX(:,k+1:end);
            gnProb.x_0             = snResult.multiMin.xOpt(:,1);
	    gnProb.LOCAL           = 1;
	 elseif dXsny < 1E-10 
            % Singular point, try disturbing x
            ss                   = ones(d,1);
            ss(rand(d,1) <= 0.5) = -1;
            h                    = 1E-6;
            %xprint(min_sn_y,'x:')
            x                    = min_sn_y + h.*ss;
            fx                   = nlp_f(x, gnProb, varargin{:});
            %fprintf('rbfSolve:   f(x) singular, disturb f(x) = %e\n',fx);
            gnProb.x_0           = x;
         else
            gnProb.x_0       = min_sn_y;
         end
         if size(snResult.multiMin.xOpt,2) == 1
            gnProb.X0        = [];
	 else
            gnProb.X0        = snResult.multiMin.xOpt(:,2:end);
         end
         gnProb.GO.ProbL.CGO = gnProb.CGO;
         gnResult            = gnSolve(gnProb);
         xInf                = gnResult.xLoc;
         fInf                = gnResult.fLoc;
         gnProb.CGO.modN     = modN;
         gnProb.CGO.fnStar   = fnStar;
         gnProb.CGO.alpha    = alpha;
         NLP_x=[]; NLP_f=[]; NARG = [];
      end
      if PriLev > 3
         fprintf('Local search OK, no need switch to target value minimization of gn_f\n')
         fprintf('min ||min_sn_y-X|| = %f\n',dXsny)
%HKHUS - Here a local integer procedure is needed
      end
      if AddSurfMin > 0
         % Pick up global multiMin solutions for analysis of option AddSurfMin
         IX   = snResult.multiMin.IX;
         xOpt = snResult.multiMin.xOpt;
         fOpt = snResult.multiMin.fOpt;
      end
      % HKH - NOTE, hard wired 3 ways of doing last cycle step
      % OLDuc = 1  Use global optimization (normal use)
      % OLDuc = 2  Do local optimization with GOlocalSolver, call snSolve
      % OLDuc = 3  Do local optimization with localSolver, call tomRun
      OLDuc = 1;
      if OLDuc == 1
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
         if OLDuc == 2
            snProb.CGO.globalSolver = GOlocalSolver;
            snProb.CGO.localSolver  = GOlocalSolver;
            snProb.CGO.fnStar       = NaN;
            snProb.x_0              = x_min;
            snProb.PriLev           = PriSub;
            snProb.GO.ProbL.CGO     = snProb.CGO; % For safety, snSolve calls local with snProb
            snResult                = snSolve(snProb);
         else
            ProbL                   = snProb.GO.ProbL;
            ProbL.CGO.fnStar        = NaN;
            ProbL.x_0               = x_min;
            if ~isempty(IntVars)
               x_LLx = snProb.x_L(IntVars);
               x_UUx = snProb.x_U(IntVars);
            end
            if ~isempty(IntVars)
               % Fix integer variables values to nearest integer point
               x00 = max(x_LLx,min(x_UUx,round(ProbL.x_0(IntVars))));
               ProbL.x_0(IntVars) = x00;
               ProbL.x_L(IntVars) = x00;
               ProbL.x_U(IntVars) = x00;
            end
            snResult                = tomRun(localSolver,ProbL,2);
            if ~isempty(IntVars)
               % Reset bounds
               ProbL.x_L(IntVars) = x_LLx;
               ProbL.x_U(IntVars) = x_UUx;
            end
         end
         fLoc                    = snResult.f_k;
         xLoc                    = snResult.x_k(:,1);
         dXsny                   = min(tomsol(30,xLoc,X));
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
                     fprintf('LipUpp %8.2g ',LipUpp);
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
                  if LipfM <= LipUpp & distIX > 0.01
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
                  xprinti(ix)
                  disp(xOL)
               end
            end
            addlocStep = length(fOL);
         end
      end
      if dXsny == 0, xLoc = addCloseX(xLoc,X,x_LL,x_UU,IntVars); end
   end
   if modN >= -1 & (modN ~= N | (modN == N & ~ucOK))

      if modN == -1 | SolveInf
         NLP_x=[]; NLP_f=[]; NARG = [];
         gnProb.CGO.fnStar   = 0;
         gnProb.CGO.alpha    = 0;
         gnProb.CGO.modN     = -1;
         gnProb.CGO.fnStar   = 0;
         gnProb.CGO.alpha    = 0;
         gnProb.IX           = [];
         gnProb.xOpt         = [];
         gnProb.fOpt         = [];
         gnProb.PriLev       = PriSub;
         dXsny               = min(tomsol(30,min_sn_y,X));
         if dXsny < 1E-10 & length(IntVars) == d
            k  = size(gnProb.CGO.X,2);
            XX = mipFeasible(gnProb, 4000, 1, 1, 0);
	    snResult.multiMin.xOpt = XX(:,k+1:end);
            gnProb.x_0             = snResult.multiMin.xOpt(:,1);
	    gnProb.LOCAL           = 1;
	 elseif dXsny < 1E-10 
            % Singular point, try disturbing x
            ss                   = ones(d,1);
            ss(rand(d,1) <= 0.5) = -1;
            h                    = 1E-6;
            %xprint(min_sn_y,'x:')
            x                    = min_sn_y + h.*ss;
            fx                   = nlp_f(x, gnProb, varargin{:});
            %fprintf('rbfSolve:   f(x) singular, disturb f(x) = %e\n',fx);
            gnProb.x_0           = x;
         else
            gnProb.x_0       = min_sn_y;
         end
         if size(snResult.multiMin.xOpt,2) == 1
            gnProb.X0        = [];
	 else
            gnProb.X0        = snResult.multiMin.xOpt(:,2:end);
         end
         gnProb.GO.ProbL.CGO = gnProb.CGO;
         gnResult            = gnSolve(gnProb);
         xLoc                = gnResult.xLoc;
         fLoc                = gnResult.fLoc;
         xInf                = xLoc;
         fInf                = fLoc;
         gnProb.CGO.modN     = modN;
         gnProb.CGO.fnStar   = fnStar;
         gnProb.CGO.alpha    = alpha;
         NLP_x=[]; NLP_f=[]; NARG = [];
      end
      if modN ~= -1
         gnProb.CGO.modN     = modN;
         gnProb.PriLev       = PriSub;
         if fnStar > fMin & PriLev > 1
            fprintf('------ WARNING!!! fnStar > fMin)');
            fprintf('fnStar %18.12f ',fnStar);
            fprintf('fMin %18.12f ',fMin);
            fprintf('\n');
         end

         % Set iteration dependent parameters used in gnSolve, output from snSolve
         if TargetMin == 2
            gnProb.IX   = snResult.multiMin.IX;
            gnProb.xOpt = snResult.multiMin.xOpt;
            gnProb.fOpt = snResult.multiMin.fOpt;
         end
         dXsny               = min(tomsol(30,min_sn_y,X));
         if dXsny < 1E-10 & length(IntVars) == d
            k  = size(gnProb.CGO.X,2);
            XX = mipFeasible(gnProb, 4000, 1, 1, 0);
	    snResult.multiMin.xOpt = XX(:,k+1:end);
            gnProb.x_0             = snResult.multiMin.xOpt(:,1);
	    % gnProb.LOCAL           = 1;
	 elseif dXsny < 1E-10 
            % Singular point, try disturbing x
            ss                   = ones(d,1);
            ss(rand(d,1) <= 0.5) = -1;
            h                    = 1E-6;
            %xprint(min_sn_y,'x:')
            x                    = min_sn_y + h.*ss;
            fx                   = nlp_f(x, gnProb, varargin{:});
            %fprintf('rbfSolve:   f(x) singular, disturb f(x) = %e\n',fx);
            gnProb.x_0           = x;
         else
            gnProb.x_0       = min_sn_y;
         end
         if size(snResult.multiMin.xOpt,2) == 1
            gnProb.X0        = [];
	 else
            gnProb.X0        = snResult.multiMin.xOpt(:,2:end);
         end
         gnProb.GO.ProbL.CGO = gnProb.CGO;
         gnResult            = gnSolve(gnProb);
         xLoc                = gnResult.xLoc;
         fLoc                = gnResult.fLoc;
      end
   end

   % Best point found on surface is xLoc - set as xNew
   xNew = xLoc;
   ix = find(all(X==xNew*ones(1,size(X,2))));
   if ~isempty(ix) & ~ucOK & modN == N
      xLoc = addCloseX(xLoc,X,x_LL,x_UU,IntVars);
      xNew = xLoc;
      ix = find(all(X==xNew*ones(1,size(X,2))));
   end
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
         fPen = fNew;
         if dLin > 0
            L = Prob.A*O_new;
            fPen = fPen+sum(max(0,max(Prob.b_L-bTol*max(1,abs(Prob.b_L))-L,...
	                              L-bTol*max(1,abs(Prob.b_U))-Prob.b_U)));
         end
         if dCon(1) > 0 
            C = cgo_c(O_new, Prob, varargin{:});
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

     % Are we trying to remove noise?
     if(NOISE == 1)
        % Update noise column
        updateNoiseColumn(etas, Prob);
     end
     
     cgolib(112, Prob.CGOLIB.obs, Prob.CGOLIB.rbfcol);
     
     % The surfaces need to be interpolated again.
     % NOTE: If not all surfaces are used in all iterations, then it could
     %       be better not to interpolate all surfaces here, but only
     %       interpolate those that will be used soon.
     
     if 1
         mindist = cgolib(114, Prob.CGOLIB.obs);
%         for i = 1:length(Prob.CGOLIB.rbftype)
          for i = rbfType
          % Set gamma for the gamma-dependent surfaces
             if(i == 3 || i == 4)
                 options = struct('');
                 options(1).PHIGAMMA = sqrt(mindist);
                 cgolib(350, Prob.CGOLIB.rbftype(i), options);
             elseif(i == 5)
                 options = struct('');
                 options(1).PHIGAMMA = 1/mindist;
                 cgolib(350, Prob.CGOLIB.rbftype(i), options);
             end
             
             % Interpolate
             cgolib(302, Prob.CGOLIB.rbftype(i));
             if(cgolib(323, Prob.CGOLIB.rbftype(i)) ~= 0)
                 % The matrix is probably very ill-conditioned.    
             end
         end
     else
         mindist = cgolib(114, Prob.CGOLIB.obs);
         fprintf('------------------- NY PUNKT -------------------\n');
         disp(['Minsta avstånd är nu: ' num2str(mindist)]);
         fprintf('\n');
         for i = 1:length(Prob.CGOLIB.rbftype)
             if(i == 3 || i == 4)
                 options = struct('');
                 options(1).PHIGAMMA = sqrt(mindist);
                 cgolib(350, Prob.CGOLIB.rbftype(i), options);
             end
             if(i == 5)
                 options = struct('');
                 options(1).PHIGAMMA = 1/mindist;
                 cgolib(350, Prob.CGOLIB.rbftype(i), options); 
             end
             fprintf('RBF typ: %i\n', i);

             cgolib(302, Prob.CGOLIB.rbftype(i));
             fprintf('.... kond_10: %f\n', log(cond(cgolib(325,Prob.CGOLIB.rbftype(i))))/log(10));

             if(cgolib(323, Prob.CGOLIB.rbftype(i)) ~= 0)
                 % This seems to occur for surface type 3, 4 and 5. This must be
                 % investigated.
                 %fprintf(['Warning: Interpolation of rbf surface type ' num2str(i) ' might be bad.\n']);
                 hohoh=3;    
             end
             cval = cgolib(324, Prob.CGOLIB.rbftype(i));
             fprintf('.... crosval: %f\n', cval);
         end
     end
    
     
     % For now, assume everything was ok
     Update = 1;
     
     % F_m0 should be set if REPLACE == 1. 
     if(REPLACE >= 1)
       nrows = cgolib(107, Prob.CGOLIB.obs);
       
       % FRHE: Problem here: If we have noise, F_m0 changes in each iteration. 
       % Should F_m be updated even for other REPLACE when there's noise?
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
      Fpen  = [Fpen;fPen];
      F00   = [F00;f00New];
      if ~isempty(CcNew)
         Cc    = [Cc,CcNew];
      end
      X     = [X xNew];
      O     = [O O_new];
      O_pre = O_new;
      n     = n + 1;

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
         FLOW     = nFunc;
      end
      % fprintf('Diff f and fPen %20.16f\n',fNew-fPen)
      NOUPDATE = 0;
   elseif Update == -2
      % Update == -2 New point bad, even refactorization failed
      % Infeasible problem, ill-conditioning?

      Dist = tomsol(30,X(:,fIdx),X);
      Dist(fIdx) = Inf;

      control = -1;
      VALUE = 1;
      while control < 0
	if(usecgolib == 1)
	  error('This can''t happen for cgolib at the moment.');
	else
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
	  disp('make init again');
	  control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
	  VALUE = 0;
	end
      end
      if control < 0
         fprintf('New initial interpolation failed');
         tomsol(25) % Deallocates memory
         Result.ExitFlag  = 1;
         Result.ExitText  = 'New Initial interpolation failed';
         Result.CGO.fGoal = fGoal;
         Result.CGO.Its   = Its;
         Result.DIGIT     = TESTDIG;
         Result           = endSolve(Prob,Result);
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

   if isempty(fGoal) | isinf(fGoal)
      RelErr = NaN;
   else
      if fGoal == 0
         RelErr = fMin-fGoal;
      else
         RelErr = (fMin-fGoal)/abs(fGoal);
      end
   end
   surfErr = fNew-snNew;

   Its.Iter(Iter)             = Iter;
   Its.n(Iter)                = n;
   Its.nFunc(Iter)            = nFunc;
   Its.modN(Iter)             = modN;
   if idea == 1
      Its.fnStar(Iter)        = fnStar;
   else
      Its.fnStar(Iter)        = alpha;
   end
   Its.fDiff(Iter)            = fDiff;
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
      snOptH             = nlp_H(xOptS,snProb);
      snOptE             = sum(eig(snOptH)<0);
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
      saveCGO(1,'rbfSolve',Name,O,F,X,F_m,F00,Cc,nInit,Fpen,fMinIdx,rngState);
   end
   

   if PriLev > 1 | IterPrint
      fprintf('\n');
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
         fprintf('  RelErr %10.6f', RelErr);
      end
      fprintf('\n');
      fprintf('  ');
      if ~isempty(CcMin)
         fprintf('CcMin:  ');
         fprintf('%11.8f ', CcMin);
      end
      xprint(O_new,'xNew:',' %12.8f',8)
      fprintf('  snErr %7.4f', surfErr);
      %fprintf(' min_sn %7.4f', min_sn);
      if fLoc > 1E30
         fprintf(' fLoc %7.4f', inf);
      else
         fprintf(' fLoc %7.4f', fLoc);
      end
      %NEWHKHfprintf(' %11.8f', fGlob);
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
      else
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
         fprintf('fRed %f. ',fRed)
         if modN ~= N % | (Rescue == 0 & ~ucOK)
            %fprintf('hn %f ',-1/gnResult.f_k);
            fprintf('hn %f ',-1/gnResult.f_k);
            if idea == 1
               %fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+1/gnResult.f_k);
               fprintf('HNerr %e ',myNew*(snNew-fnStar)^2+1/gnResult.f_k);
            end
         end
         fprintf('LipU %f. ',LipUpp)
         fprintf('LipL %f. ',LipLow)
         fprintf('\n');
      end
      if ~isempty(xOptS) & PriLev > 3
         fprintf('  dXO (min||xOpt-X||) %f ',dXO);
         fprintf('sn_f(xOpt) %8.5f ',snOptf);
         fprintf('||sn_g(xOpt)|| %f ',snOptg);
         fprintf('negeig(sn_H(xOpt)) %d ',snOptE);
         fprintf('\n');
      end
   end

   % ********** CONVERGENCE TEST **********
   if n >= nMax, convflag = 7; end
   if convflag == 0
      convflag = CGOisClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
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

   if abs(Update) == 1
      fDiffOld = fDiff;
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
               if dCon(1) > 0 & dCon(2) == 0
                  C = nlp_c(xF, Prob, varargin{:});
                  nCon = nCon + 1;
                  cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
               elseif dCon(1) > 0 & dCon(2) > 0
                  C = nlp_c(xF, snProb, varargin{:});
                  nCon = nCon + 1;
                  cVal = cVal+sum(max(...
		           0,max(Prob.CGO.c_L-cTol-C,C-cTol-Prob.CGO.c_U)));
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
                  if dCon(1) > 0 & dCon(2) == 0
                     C = nlp_c(xF, Prob, varargin{:});
                     nCon = nCon + 1;
                     cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
                   elseif dCon(1) > 0 & dCon(2) > 0
                     C = nlp_c(xF, snProb, varargin{:});
                     nCon = nCon + 1;
                     cVal = cVal+sum(max(...
		            0,max(Prob.CGO.c_L-cTol-C,C-cTol-Prob.CGO.c_U)));
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
                     if dCon(1) > 0 & dCon(2) == 0
                        C = nlp_c(xF, Prob, varargin{:});
                        nCon = nCon + 1;
                        cVal = cVal+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
                     elseif dCon(1) > 0 & dCon(2) > 0
                        C = nlp_c(xF, snProb, varargin{:});
                        nCon = nCon + 1;
                        cVal = cVal+sum(max(...
		               0,max(Prob.CGO.c_L-cTol-C,C-cTol-Prob.CGO.c_U)));
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
Result.CGO.WarmStartInfo = saveCGO(2,'rbfSolve', Name,O,F,X,F_m,F00,Cc,...
                          nInit,Fpen,fMinIdx, rngState);

% All points i with F(i)=f_min
if SCALE
   Result.x_k      = tomsol(9,x_L,X(:,F==fMin),x_D);
else
   Result.x_k      = X(:,F==fMin);
end
if isempty(Result.x_k)
   if SCALE
      Result.x_k      = tomsol(9,x_L,X(:,Fpen==fMin),x_D);
   else
      Result.x_k      = X(:,Fpen==fMin);
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

if dLin > 0
   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at all x_k
end

%if dLin > 0 & size(Result.x_k,2)==1
%   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at best x_k
%end

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

% -------------- End of main rbfSolve routine -------------------
% --------------------------------------------------------------
% --------------------------------------------------------------

% --------------------------------------------------------------
function [idea, rbfType, N, infStep, fStarRule, ...
   DeltaRule, MaxCycle, eps_sn, TargetMin, AddSurfMin, ...
   usecgolib, NOISE] = getRBFProbVars(Prob, IntVars)
% --------------------------------------------------------------


if isempty(Prob.CGO)
   idea         = []; rbfType      = []; N         = []; infStep    = [];
   fStarRule    = []; DeltaRule    = []; MaxCycle  = []; eps_sn     = [];
   TargetMin    = []; AddSurfMin   = [];
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
   if isfield(Prob.CGO, 'NOISE')
     NOISE = Prob.CGO.NOISE;
   else
     NOISE = 0;
   end
end

if isempty(idea),         idea = 1; end
if isempty(rbfType),      rbfType = 2; end
if isempty(infStep),      infStep = 0; end
if isempty(fStarRule),    fStarRule = 1; end
if isempty(DeltaRule),    DeltaRule = 1; end
if isempty(MaxCycle)     
   if isempty(IntVars)
      MaxCycle = 10; 
   else
      MaxCycle = 1000; 
   end
end
if isempty(eps_sn),       eps_sn    = 1E-7; end
if isempty(AddSurfMin),   AddSurfMin = 0; end
if isempty(TargetMin),    TargetMin = 3; end

if idea == 1      % Cycle length for direct fnStar strategy
   if fStarRule == 3
      infStep = 1; % Must have infStep in strategy 3, (Gutmann I)
      if isempty(N), N = 1; end
   else
      if isempty(N), N = 4; end
   end
elseif idea == 2  % Cycle length for indirect alpha strategy
   usecgolib = 0; % gn_g does not work with alpha strategy
   %if isempty(N), N = 3; end
   N = 3; % Always N = 3
end
% --------------------------------------------------------------
function [alpha, ucOK] = getAlpha(modN, alpha, F, F_m, min_sn, fMin, ...
          f_Low, eps_sn, Update, Iter, fMinIter, fDiffOld, OldAlpha)
% --------------------------------------------------------------
ucOK  = 1;
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
         end
      case 2
         if fMinIter == Iter-1
            %alpha = min(fDiff/2,1);
            % test, working?
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
   %if PriLev > -1
       fprintf('Alpha same as previous step, adjust to %10.3f\n',alpha);
   %end
end

% Comments on NOISE
%
% This is what's been done:
%  A data column (Prob.CGOLIB.noisecol) is created.
%  An RBF surface (Prob.CGOLIB.noiserbf) is created based on that column.
%  
%  A function, updateNoiseColumn, has been written. Given a set
%  of etas, it computes the error for each eta and makes a cross-
%  validation to get the eta that gives the smallest error. The
%  corrected funtion values (f+e) are put into noisecol.
%  
%  It makes this update every iteration.
%
%  If Prob.CGO.NOISE is set to 1, then the REPLACE transform
%  columns use the noisecol as source instead of datacol. So,
%  the noise is corrected first, THEN the replaces and transoforms
%  are applied.	
%
% Issues:
%
%  No correction is done to Fm, except in the case when REPLACE >=1.
%  (Because in that case, corrections were done before)
%  Fm is changing in every iteration within CGOLIB, as the replaces
%  and transforms are functions of the noise corrected column. This
%  column is changed every iteration, so Fm(i) doesn't stay the
%  same between iterations. So: The F-vectors must be updated 
%  correctly. An error that may occur because of this is that 
%  min_sn > fMinF and the solver retries with another solver.
%
%  eta seems to stick with high values. When eta = 1, the rbf surface
%  is simply a first degree polynomial and probably doesn't predict
%  the function good.
%
%  Better methods for finding eta would be good. Now there is a
%  list of etas, and each of them is tried, and the best one is
%  taken. It is an optimization problem, and the objective is costly,
%  as it is a cross-validation computation.
%
%  Now the default RBF type is used for computing the errors (don't
%  remember if it is cubic or thin-plate). Maybe the RBF type of
%  noiserbf should be the same as the one the user chose!
%  
%  I don't know if this is an NOISE issue, but for some problems,
%  points are sampled too close to each other. This has an impact
%  on the NOISE-mechanisn, as we want to compute the inverse of A.
%
%  In tests, the NOISE-mechanism doesn't seem to work very well.
%  It usually works better without it. Maybe something is wrong
%  in the implementation somewhere. But I also have a theory why
%  it doens't prove well in the tests:
%   What the NOISE-mechanism does is to smooth out the function.
%   Sharp minimums may disappear in the smoothing. Maybe this
%   method works better for flat minimums? 

function updateNoiseColumn(etas, Prob)

rbf = Prob.CGOLIB.noiserbf;
n = cgolib(107, Prob.CGOLIB.obs);

cgolib(302,rbf);
A = cgolib(325,rbf);
Ainv = inv(A);
B = Ainv(1:n, 1:n);

mindist = cgolib(114, Prob.CGOLIB.obs);
if(mindist < 1e-7)
    fprintf('UPDATENOISECOLUMN: Minimum distance = %e\n', mindist);
end

BAB = B'*A(1:n,1:n)*B;

% Get original f-values
f = cgolib(110, Prob.CGOLIB.obs, n, 1, [], Prob.CGOLIB.datacol);
rhs = BAB*f;

beste = zeros(n,1);
bestdelta = 1e300;
delta = 1e300;
besteta = etas(1);

% Save the deltas. One can plot them. In all cases
% I have been looking, delta has had only one minimum.
deltai = [];
for i = 1:length(etas)
   
    eta = etas(i);

    C = ((eta-1)/eta*eye(n)-BAB);

    % Solve for the errors
    e = C\rhs;

    % Set the noise column to f+e
    cgolib(106, Prob.CGOLIB.obs, [], Prob.CGOLIB.noisecol, f+e);

    if(length(etas) > 1)
        % Interpolte surface
        cgolib(302, Prob.CGOLIB.noiserbf);
    
        % Cross-validate
        delta = cgolib(324, Prob.CGOLIB.noiserbf, f);
        if(delta < bestdelta)
            beste = e;
            bestdelta = delta;
            besteta = eta;
        end
        deltai(i) = delta;
    end
end

if(delta ~= bestdelta)
    cgolib(106, Prob.CGOLIB.obs, [], Prob.CGOLIB.noisecol, f+beste);
end
    
fprintf('Best eta: %f\n', besteta);    

% % % % --------------------------------------------------------------
% % % function PLOTCRAP(Prob, x_L, x_U)
% % % % --------------------------------------------------------------
% % % % Just save code for plotting, for now
% % %    % ********** PLOTTING **********
% % %    if PLOT & d > 1
% % %
% % %       if 10 & ~(modN == N)
% % %          switch lower(globalSolver)
% % %          case 'glcdirect'
% % %               glb=load('glcDirectSave.mat','C');
% % %          case 'glbdirect'
% % %               glb=load('glbDirectSave.mat','C');
% % %          case 'glcfast'
% % %               glb=load('glcFastSave.mat','C');
% % %          case 'glbfast'
% % %               glb=load('glbFastSave.mat','C');
% % %          case 'glbsolve'
% % %               glb=load('glbSave.mat','C');
% % %          case 'glcsolve'
% % %               glb=load('glcSave.mat','C');
% % %          otherwise
% % %               %break;
% % %          end
% % %          % Transform to original space
% % %          if SCALE
% % %             CCC   = tomsol(9, x_L, glb.C, x_D);
% % %          else
% % %             CCC   = glb.C;
% % %          end
% % %          CCC = zeros(size(glb.C));
% % %
% % %          plot(CCC(1,:),CCC(2,:),'.r');
% % %          hold on
% % %          plot(O(1,:),O(2,:),'*');
% % %       else
% % %          plot(O(1,:),O(2,:),'*');
% % %          hold on
% % %          TT = 0:0.1:1;
% % %          YY(1,:) = c0(1)+TT*(O_new(1)-c0(1));
% % %          YY(2,:) = c0(2)+TT*(O_new(2)-c0(2));
% % %          plot(YY(1,:),YY(2,:),'-g')
% % %       end
% % %       plot(O_new(1),O_new(2),'*g');
% % %       x_opt = Prob.x_opt;
% % %       if ~isempty(x_opt)
% % %          if min(size(x_opt))==1
% % %             x_opt = x_opt(:);
% % %          end
% % %          plot(x_opt(1,:),x_opt(2,:),'*k');
% % %       end
% % %       hold off
% % %       grid on
% % %       zoom on
% % %    end % if plot
% % %    if PLOT & d == 1
% % %          switch lower(globalSolver)
% % %          case 'glcdirect'
% % %               glb=load('glcDirectSave.mat','C');
% % %          case 'glbdirect'
% % %               glb=load('glbDirectSave.mat','C');
% % %          case 'glcfast'
% % %               glb=load('glcFastSave.mat','C');
% % %          case 'glbfast'
% % %               glb=load('glbFastSave.mat','C');
% % %          case 'glbsolve'
% % %               glb=load('glbSave.mat','C');
% % %          case 'glcsolve'
% % %               glb=load('glcSave.mat','C');
% % %          otherwise
% % %               %break;
% % %          end
% % %          % Transform to original space
% % %          if SCALE
% % %             CCC   = tomsol(9, x_L, glb.C, x_D);
% % %          else
% % %             CCC   = glb.C;
% % %          end
% % %          CCC = zeros(size(glb.C));
% % %       plot(CCC(1,:),zeros(size(CCC,2)),'.r');
% % %       plot(O_new(1),fLoc,'*g');
% % %       x_opt = Prob.x_opt;
% % %       if ~isempty(x_opt)
% % %          if min(size(x_opt))==1
% % %             x_opt = x_opt(:);
% % %          end
% % %          plot(x_opt(1,:),x_opt(2,:),'*k');
% % %       end
% % %       hold off
% % %       grid on
% % %       zoom on
% % %    end


% MODIFICATION LOG
%
% 000112  mbk  First version written
% 001101  hkh  Revision for new tests, speedups
% 010431  hkh  Total revision, making it similar to glbSolve
% 010723  hkh  Send MaxFunc to Mex
% 011110  hkh  Fixes for GUI. Revision.
% 011206  hkh  Alternative cycle strategies are tried.
% 020103  hkh  Use general field CGO and cgoSave.mat for warm start
% 020614  hkh  Correcting errors in print statements
% 020820  hkh  Print sub problem results dependent on PriLev
% 020820  hkh  Try local search for all possible global solutions
% 020823  hkh  Error in transforming linear constraints if SCALE = 1
% 020824  hkh  Accept local solution only if linear constraints fulfilled
% 020828  hkh  Use GetSolver to select default local solver
% 020831  hkh  Add parameters in Prob.GO as input parameters
% 020831  hkh  Add CGO.LOCAL, cycle length CGO.N as input parameters
% 020831  hkh  Add CGO.infStep, add inf target value step in cycle
% 020831  hkh  Add CGO.AddMP, add midpoint to corner strategy if true
% 020831  hkh  Add CGO.fStarRule and CGO.DeltaRule
% 020831  hkh  Added gutmann d+1 initial value strategy
% 020901  hkh  Create two more fStar search strategies, now all Gutmann
% 020904  hkh  Set GOMaxFunc and GOMaxIter empty if not initialized.
% 021020  hkh  MaxFunc and MaxIter wrongly set empty if Prob.GO empty
% 021020  hkh  MaxFunc and MaxIter wrongly set empty if Prob.GO empty
% 021020  hkh  Use Fpen=F+sum(max(0,dev)) to handle infeasible initial points
% 030221  hkh  isempty(GOLocalSolver) should be isempty(GO.LocalSolver)
% 030222  hkh  Add check for empty Result.x_k, use best infeasible solution
% 030223  hkh  randomtest used size as variable name, changed to m
% 030223  hkh  Changed variable names, (xNew,fNew) for f(x) point
%              (xGlob,fGlob) for all global opt solutions on surface,
%              (xLoc,fLoc) for best solution found on surface
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin in snProb and gnProb
% 040226  hkh  Set field Prob.MIP in gnProb.MIP and snProb.MIP
% 040226  hkh  Set globalSolver='glcCluster' if Prob.MIP.IntVars is non-empty
% 040303  hkh  Making rbfSolve handle mixed-integer problems
% 040306  hkh  Add CGO.RandState, to initialize the random generator
% 040306  hkh  Check for cycle of same points
% 040306  hkh  Set Result.Inform to convflag
% 040307  hkh  Revised functions gutmann and corners, added new sampleInts
% 040307  hkh  Revise feasibility handling, print fMinI, fMinF, use Feasible
% 040308  hkh  Special treatment of pure IP, with glcSolve and glcFast
% 040318  hkh  Major tuning of algorithm, Update codes etc.
% 040318  hkh  Set CGO.F as column always. Revised handling of initial F,X
% 040319  hkh  New optional input CGO.CX for initial constraint values
% 040321  hkh  Set initial FLOWER to fIdx(1), not nFunc
% 040323  hkh  Set cTol and bTol in suboptimizations, avoiding infeasibility
% 040323  hkh  Default for Percent if dLin+dCon > 0, use -5000 instead of -d
% 040404  hkh  Major revision, global solvers do not return empty
% 040406  hkh  Save fMinIdx (=fIdx(1)) on cgoSave.mat
% 040412  hkh  Add calls to ProbCheck
% 040414  hkh  Do not set MENU or LargeScale in gnProb and snProb
% 040414  hkh  Return Result.c_k and Result.Ax
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 040928  hkh  sampleInts, x_L(j), offset from 0, must be added
% 040928  hkh  If localSolver is glcCluster, GOMaxFunc must be used
% 041005  hkh  Use abs(fMin-min_sn) <= eps_sn*abs(fMin), not fMin-min_sn
% 041005  hkh  Revise STILL TOO CLOSE handling if glcCluster
% 041005  hkh  Adding ExitFlag=8 for no progress for MaxCycle*(N+1) func evals
% 041005  hkh  CGO.MaxCycle new input, default 10
% 041005  hkh  Add localSolver and globalSolver to SolverAlgorithm output
% 041006  hkh  If Unc Sol too close, use global solver on snProb, not gnProb
% 041018  hkh  Set snProb.GO = Prob.GO, because global solver sometimes used
% 041019  hkh  Avoid global step if "TOO CLOSE", set eps_sn = 1E-6
% 041114  hkh  Add parameters for solver minlpSolve
% 041114  hkh  Rerun with rounded IntVars if local solution not int feasible
% 041114  hkh  Change default strategy to idea 1
% 041123  hkh  Change call to tomRun
% 041130  hkh  Remove unused code (previously used when ~ucOK and modN==N)
% 050216  hkh  New function daceInts, round experimental design points
% 050216  hkh  Add check on duplicate points in sampleInts and daceInts
% 050216  hkh  Use daceInts instead of sampleInts when for DACE strategies
% 050216  hkh  Add new output field Result.CGO with tuning params used
% 050321  hkh  Forgot to define n = size(X,2) when changing initial strategy
% 050322  hkh  Pure IP: Test if all integer values tried (nMax == n)
% 050322  hkh  Wrong check for duplicates in daceInts, also round only IntVars
% 050322  hkh  Test if all integer values tried, if all continuous fixed
% 050322  hkh  Add test if all variables are fixed
% 050423  hkh  Set field gnProb.FUNCS.d2c to rbf_d2c
% 050426  hkh  If initial interpolation fails & REPLACE=0, retry with REPLACE=1
% 050426  hkh  Print median(F),min(F),comments if initial interpolation fails
% 050426  hkh  Improved ellipsoid initial method. Use tomsol distance
% 050506  hkh  Correction of infStep method
% 051012  hkh  Wrong call to rand if RandState < 0
% 051012  hkh  Avoid estimation of gradient in endSolve setting GradEv = -1
% 051015  hkh  max_F computed wrong, from F, not F_sort if REPLACE=0
% 060126  hkh  Changed glcFast to glcDirect as default, added glcDirect setup
% 060206  hkh  Add new input nSample and eps_sn
% 060207  hkh  Use new function expDesign for initial experimental design
% 060508  hkh  Add output field Result.Dig3 for Iter to 3 digits convergence
% 060508  hkh  Add full comments on IterPrint printout
% 060508  hkh  If xNew in X, use min_sn_y as RESCUE, print R after Cycle
% 060626  hkh  Check if Prob.x_opt empty before access
% 060626  hkh  Make a structure of Result.Dig3 with fields Iter,FuncEv,CPUtime
% 060629  hkh  Add DIRECT methods as initial design methods, Percent > 100
% 060629  hkh  Changed default for eps_sn to 1E-7
% 060814  med  FUNCS used for callbacks instead
% 060816  hkh  Prob.MIP.IntVars also set if isempty(Prob.MIP), expDesign crash
% 060817  hkh  Test if finite bounds. Define Result before tests
% 061124  hkh  Lower letter in case tests. Add multiMin support
% 061127  hkh  Default cycle length N now 5 not 4, Gutmann used 2,3,4 only
% 061127  hkh  Add output field Result.DIGIT for 1,...,6 digits convergence
% 070221  hkh  Set upper bound on MaxFunc to 5000, avoid crash in mex
% 070222  hkh  Revise IntVars handling, use new format
% 070223  hkh  Safe guard REPLACE, other solvers use more than 0,1
% 070301  hkh  Define some fields in Result structure when interpolation failure
% 070308  hkh  Define backupSolver1 and backupSolver2, safe snSolve computation
% 070925  hkh  Set snProb.LOCAL=LOCAL, was set in snSolve in wrong way
% 071002  hkh  Now 3 exp-designs based on adjacent corners,-997,-998,-999
% 071002  hkh  Add new feature, and input parameter AddSurfMin
% 071006  hkh  Do init of random generator only if no warm start, NaN no init
% 071007  hkh  Save random generator state in rngState, reinit rand if WarmStart
% 071007  hkh  Add TargetMin parameter, 4 strategies to select minimum
% 071007  hkh  Add f_Low in Result.CGO, has influence to fDiff if > -inf
% 071007  hkh  Input Prob.PriLevSub for output in subproblem solvers snSolve,gnSolve
% 071010  hkh  Use gnProb.x_0 and gnProb.X0 as additional initial points
% 080410  hkh  Complete revision in modules, new experimental design
% 080415  hkh  Prob.GO input in Result.CGO, change CGOGLobalProb call
% 080415  hkh  Add nTrial, CLHMethod to Result.CGO
% 080416  hkh  Avoid 2 identical initpoints in local search, split in X0 and x_0
% 080417  hkh  Use relative constraint violation in fPen computation
% 080616  frhe CGOLIB calls added
% 080617  frhe Made REPLACE>1 transformation monotonic
% 080617  frhe log10 used in REPLACE>1 in accordance with help
% 080625  hkh  Check if min_sn_y is very close to sampled point, disturb
% 080626  hkh  Use Prob.PriLevSub in snSolve/gnSolve
% 080630  hkh  Use cgo_fc instead of nlp_f for costly f(x) (and costly c(x))
% 080702  frhe All possible RBF types maintained during run for demo
% 080711  hkh  Use routine saveCGO for warm start info saving for all solvers
% 080711  hkh  Save warm start info every iteration in cgoSave.mat (cgoSave1)
% 080711  hkh  Avoid estimation of Hessian in endSolve setting HessEv = -1
% 080711  hkh  Save warm start info/cgoSave if initial interpolation fails
% 080801  frhe NOISE added.
% 081105  hkh  Make getProbVars solver dependent
% 081105  hkh  Add input/output nTrial in call to getProbVars
% 090429  hkh  MaxCycle = 1000 if any IntVars, otherwise default 10 as before
% 090501  hkh  Only compute 1 set of RBF interpolation surfaces
% 090813  med  repmat call corrected
% 090824  hkh  Minor mlint revision
% 091006  hkh  Write internal function getAlpha to compute alpha in cycle
% 091006  hkh  Write external function getfnStar to compute fnStar in cycle

