% ego.m
%
% function Result = ego(Prob, varargin)
%
% ego implements the algorithm EGO by D. R. Jones,
% Matthias Schonlau and William J. Welch presented in the paper
% "Efficient Global Optimization of Expensive Black-Box Functions".
%
% The algorithm is expanded by Prof. K. Holmstrom to handle linear,
% nonlinear and integer constraints
%
% ego solves problems of the form:
%
%    min   f(x/cg
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
% function Result = ego(Prob, varargin)
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say EGOF
%   FUNCS.c   The routine to compute the nonlinear constraint, say EGOC
%             A call to tomFiles.m or glcAssign.m sets these fields.
%   x_L       Lower bounds for each element in x. Must be finite.
%   x_U       Upper bounds for each element in x. Must be finite.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   WarmStart If true, >0, ego reads the output from the last run
%             from the mat-file cgoSave.mat, and continues from the last run.
%             If Prob.CGO.WarmStartInfo has been defined through a call to
%             WarmDefGLOBAL, this field is used instead of the cgoSave.mat file
%
%   user      User field used to send information to low-level functions
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%   PriLevOpt Print Level (PriLev below)
%             0 = silent. 1 = Summary 2 = Printing each iteration
%             3 = Info about local / global solution 4 = Progress in x
%
% -----------------------
% Fields in Prob.optParam
% -----------------------
%             Defines optimization parameters. Fields used:
%  MaxFunc    Maximal number of function evaluations in ego, default 200
%             If WarmStart == 1 and MaxFunc <= nFunc (Number of f(x) used)
%             then MaxFunc = MaxFunc + nFunc
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1.
%  fGoal      Goal for function value, if empty not used
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
%              ---- Special EGO algorithm parameters in Prob.CGO ----
% EGOAlg       Main algorithm in the EGO solver (default EGOAlg == 1)
%              =1 Run expected improvement steps (modN=0,1,2,...).
%              If no f(x) improvement, use DACE surface minimum  (modN=-1) in 1 step
%
%              =2 Run expected improvement steps (modN=0) until ExpI/|fMin| < TolExpI
%              for 3 successive steps (modN=1,2,3) without f(x) improvement (fRed <=0),
%              After 2 such steps (when modN=2), 1 step using the DACE surface minimum
%              (modN=-1) is tried. If then fRed >TolExpI, reset to modN=0 steps.
%
% pEst         Norm parameters, fixed or estimated, see p0, pLow, pUpp
%              0  = Fixed constant p-value for all components (Default, p0=1.99)
%              1  = Estimate one p-value valid for all components
%              >1 = Estimate d || ||_p parameters, 1 for each component
%
% p0           Fixed p-value (pEst==0, default = 1.99) or
%              initial p-value (pEst == 1, default 1.99) or
%              d-vector of initial p-values (pEst > 1, default 1.99*ones(d,1))
%
% pLow         If pEst == 0, not used
%              if pEst == 1, lower bound on p-value (default 1.0)
%              if pEst >  1, lower bounds on p (default ones(d,1))
%
% pUpp         If pEst == 0, not used
%              if pEst == 1, upper bound on p-value (default 1.9999)
%              if pEst >  1, upper bounds on p (default 1.9999*ones(d,1))
%
% TRANSFORM    Function value transformation
%              TRANSFORM = 0  % No transformation made.
%              TRANSFORM = 1  % Median value transformation. Use
%                               REPLACE instead.
%              TRANSFORM = 2  % log(y) transformation made.
%              TRANSFORM = 3  % -log(-y) transformation made.
%              TRANSFORM = 4  % -1/y transformation made.
%              Default EGO is computing the best possible transformation
%              from the initial set of data.
%              Note! No check is made on illegal y if user gives TRANSFORM
%
% EITRANSFORM  Transformation on expected improvement function.
%              TRANSFORM = 0  % No transformation made.
%              TRANSFORM = 1  % -log(-f) transformation made.
%              TRANSFORM = 2  % -1/f transformation made.
%              Default is 1.
%
% TolExpI      Convergence tolerance for expected improvement, default 1E-7
%
% SAMPLEF      Sample criterion function:
%              0 = Expected improvment (default)
%              1 = Kushner's criterion (related option: KEPS)
%              2 = Lower confidence bounding (related option: LCBB)
%              3 = Generalized expected improvement (related option: GEIG)
%              4 = Maximum variance
%              5 = Watson and Barnes 2
%
% KEPS         The epsilon parameter in the Kushner's criterion.
%              If set to a positive value, then epsilon is taken
%              as KEPS. If set to a negative value, then epsilon
%              is taken as |KEPS|*f_min.
%              Default: -0.01
%
% GEIG         The exponent g in the generalized expected improvement
%              function.
%              Default: 2
%
% LCBB         Lower Confidence Bounding parameter b.
%              Default: 2
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
% varargin    Additional parameters to ego are sent to the costly f(x)
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
%  ExitFlag Always 0
%  Inform   0  = Normal termination
%           1  = Function value f(x) is less than fGoal
%           2  = Error in function value f(x), abs(f-fGoal) <= fTol, fGoal=0
%           3  = Relative Error in function value f(x) is less than fTol, i.e.
%                abs(f-fGoal)/abs(fGoal) <= fTol
%           4  = No new point sampled for N iteration steps
%           5  = All sample points same as the best point for N last iterations
%           6  = All sample points same as previous point for N last iterations
%           7  = All feasible integers tried
%           9  = Max CPU Time reached
%           10 = Expected improvement low for three iterations
%
%
% To make a warm start possible, ego saves the following information in
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
% Step      Type of search step
% fGoal     Goal value (if set by user)
% fMin?     Best f(x) found so far. E.g. at 27/It 12 means n=27, Iter=12
%           fMinI means the best f(x) is infeasible
%           fMinF means the best f(x) is feasible (also integer feasible)
%           IT implies reduction in last step, It no reduction last step
%
% ----------------------------------------------
% Additional information on Row 1 in iteration 0
% ----------------------------------------------
% max(F)    maximum of all f(x) in the initial set of points X
% med(F)    median of all f(x) in the initial set of points X
% rng(F)    maxF-fMin, the range of f(x) values in the initial set X
% pDist     The size of the simply bounded region, ||x_U-x_L||_2
% --------------------------------------------------------------
% Additional information on Row 2 in iteration 0 (if PriLev > 2)
% --------------------------------------------------------------
% ||x_U-x_L|| (pDist) The size of the simply bounded region (if not on Row 1)
% LipU      Maximum Lipschitz constant for initial set X
% LipL      Minimum Lipschitz constant for initial set X
% ----------------------------------------------------
% Additional information on Row 1 in iteration 1,2,...
% ----------------------------------------------------
% fNew      f(xNew), the costly function value for point tried in current Iter
% RelE      Relative distance to known or assumed global minimum (Prob.x_opt)
% DErr      (Costly f(x) value - Surface value at x) / |f(x)|, i.e.
%           Relative surface error (Actual value - Predicted value)/|Actual value|
%
% -------------------------------
% --- Row 2 Special EGO printout:
% -------------------------------
% fE        Maximum Likelihood, estimating DACE fitting parameters (Step >=0)
% fE        DACE surface minimum (Step ==-1)
% fPred     Predicted new f(x), yMin-ExpI, in transformed space (Step >=0)
% I%        Optimal value of expected improvement (ExpI) (in %) (Step >=0)
% Q%        100*abs(ExpI/yMin)), if yMin ~= 0 (Step >=0)
%           The following 7 items, if PriLev > 2:
% [.]       Number of variables x active on bound at new point xNew
%           Distances from xNew to:
% doX       Minimal distance to closest point of all previous sampled points X
% doM       Distance to current best point found xMin, f(xMin) = fMinF
% doO       Distance to global optimum (if Prob.x_opt specified)
% fRed      Function reduction fMin-fNew in new step (negative fRed - increase)
% LipU      Maximum Lipschitz constant for initial set X
% LipL      Minimum Lipschitz constant for initial set X
%           if TRANSFORM ~= 0
% yMin      fMin, transformed back to original function space
% yNew      fNew, transformed back to original function space
% yPred     yMin-ExpI, in original function space
%
% =========================================
% -------------------------------
% --- Row 3,4,5 in iteration 0
% -------------------------------
% xMin      Best point xMin in initial set X, in original coordinates (if SCALE==1)
%
% xOpt      User-given global optimum Prob.x_opt (if defined)
%
% TRANSFORM   Transformation method choosen, transforms f(x) values
% EITRANSFORM Transformation method for expected improvement
% REPLACE     Replacement method for large f(x) values
% Internal TRANSFORM  Number used internally in CGOLIB (PriLev > 2)
% =========================================
% -------------------------------
% --- Row 3 in iteration 1,2,....
% -------------------------------
% xNew      xNew scaled back to original coordinates (if SCALE==1)
% ------------------------------
% --- Row 4 (if dimension d < 5)
% ------------------------------
% th,p      theta and p vector (length d + d) in DACE estimation
% ------------------------------
% --- Row 4 (if dimension d > 5)
% ------------------------------
% p         p vector in DACE estimation
% ------------------------------
% --- Row 5 (if dimension d > 5)
% ------------------------------
% theta     theta vector in DACE estimation
% -----------------------------------------------------------------
% --- Row 6 (if global minimum Prob.x_opt is given, and PriLev > 2
% -----------------------------------------------------------------
%
% -----------
% Iteration 1,2,... (if global optimum given in Prob.x_opt):
% -----------
% dXO                  Minimal distance from global optimum to closest
%                      point of all sampled points X in experimental design
% DACE(xOpt)           The DACE surface value at Prob.x_opt
% ||DACE-g(xOpt)||     The norm of the DACE surface gradient at Prob.x_opt
% Negeig(DACE-H(xOpt)) Number of negative eigenvalues of the DACE
%                      surface Hessian at Prob.x_opt
%
% ======
% USAGE:
% ======
%
% Let the name of the problem be "EGOF Test"
% The function EGOF is best written as
%     function f = EGOF(x, Prob)
% Then any information, say u and W is easily sent to EGOF (and the constraint
% function EGOC, if present) using the Prob structure. See the example below.
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%    Prob   = glcAssign('EGOF',x_L,x_U,'EGOF Test',A,b_L,b_U,'EGOC',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown on the two next lines
%    Prob.optParam.MaxFunc = 300; % Change max number of function evaluations
%
% Direct solver call:
%    Result = ego(Prob);
%    PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('ego',Prob,2);
%
% The user function EGOF is written as
%
%    function f = EGOF(x, Prob)
%    u = Prob.user.u; W = Prob.user.W;
%    f = "some function of x, u and W"
%
% It is also possible to use the function format
%    function f = EGOF(x)
% but then any additional parameters must be sent as global variables.
%
% The user function EGOC, computing the nonlinear constraints, is written as
%
%    function c = EGOC(x, Prob)
%    u = Prob.user.u; W = Prob.user.W;
%    c = "some vector function of x, V and W"
%
% Note! If EGOF has the 2nd input argument Prob, also EGOC must have that.
%
% To make a restart, just set the restart flag, and call ego once again:
%
%    Prob.WarmStart = 1;
%    Result = tomRun('ego',Prob,2);
%
% Another option for warm start is to put the restart information from the
% Result structure to the Prob structure and then call a CGO solver again:
%    Prob = WarmDefGLOBAL('ego',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('ego',Prob,2);
%
% It is also possible to run another CGO solver in the warm start:
%    Prob = WarmDefGLOBAL('ego',Prob,Result)
%    Prob.WarmStart = 1;
%    Result = tomRun('rbfSolve',Prob,2);
%
%
% IMPORTANT:
%
% Observe that when cancelling with CTRL+C during a run, some memory
% allocated by ego will not be deallocated. To deallocate, do:
%
%  >> clear cgolib

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2010 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written Sep 27, 1998.   Last modified Mars 26, 2010.

function Result = ego(Prob, varargin)

if nargin < 1
   error('ego needs one parameter structure Prob');
end
global NLP_x NLP_f NARG

DebugPriLev = 0;  % PriLev in subopt, i.e. gnProb.PriLevOpt, snProb.PriLevOpt

time                   = fix(clock);
solvType               = checkType('glc');
Prob                   = ProbCheck(Prob,'ego',solvType);
Prob                   = iniSolve(Prob,solvType,0,0);

Result                 = ResultDef(Prob);
Result.Solver          = 'ego';
Result.SolverAlgorithm = 'Efficient Global Optimization';

% Pick up bounds from Prob structure and check if OK
x_L       = Prob.x_L(:);             % Lower bounds
x_U       = Prob.x_U(:);             % Upper bounds

if isempty(x_L) | isempty(x_U)
   disp('ego requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'ego requires both lower and upper variable bounds';
   Result.DIGIT    = [];
   Result.CGO.Its  = [];
   Result          = endSolve(Prob,Result);
   return;
end
if ~(all(isfinite(x_L)) & all(isfinite(Prob.x_U)))
   disp('ego only solves box bounded problems, where both');
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
% PEN=Prob.CGO.PEN; % NOT USED NOW
SAVE = 1;
% ----------------------------------------------------------
% Set possibly redefined GOlocalSolver and DIRECT into Prob.GO
% ----------------------------------------------------------
Prob.GO.localSolver  = GOlocalSolver;
Prob.GO.DIRECT       = DIRECT;

% -------------------------------------------------
% Get special EGO variables from Prob.CGO structure
% -------------------------------------------------
[TRANSFORM, EITRANSFORM, EGOAlg, TolExpI, pEst, p0, pLow, pUpp, ...
        SAMPLEF, KEPS, LCBB, GEIG] = ...
   getEGOProbVars(Prob);
%getEGOProbVars(Prob, x_L, x_U);

AlgText = ['EGO Alg ' num2str(EGOAlg)];
if pEst == 0
   Result.SolverAlgorithm = [AlgText '. All p=' num2str(p0) ] ;
elseif pEst == 1
   Result.SolverAlgorithm = [AlgText '. One p optimized '] ;
else
   Result.SolverAlgorithm = [AlgText '. All p optimized '] ;
end

% -------------------------------------------------------------
% Get variables from Prob structure
% Redefine GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3
% and set back in Prob.GO
% -------------------------------------------------------------
[Prob ,MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, ...
   fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX, ...
   x_LL,x_UU, x_DD, pDist, dLin, dCon, ...
   Percent, nSample, nTrial, AddMP, REPLACE, ...
   backupSolver1, backupSolver2, xOptS ...
   ] = getProbVars(Prob, x_L,x_U,x_D, SCALE,globalSolver,IntVars,...
   Percent, nSample, nTrial, AddMP, REPLACE, ...
   GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3,'ego');

% HKH - right now do not allow REPLACE > 1 in EGO
REPLACE = min(1,REPLACE);

% ----------------------------------------------------------------------------------------
% HKH Currently not using f_Low in EGO
% ----------------------------------------------------------------------------------------

%
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
Result.CGO = struct(...
   'SCALE',SCALE, 'REPLACE',REPLACE, 'RandState',RandState,...
   'Percent',Percent,'nSample',nSample,'AddMP',AddMP, ...
   'nTrial',nTrial, 'CLHMethod',CLHMethod, ...
   'fTol',fTol, 'f_Low',f_Low,'SMOOTH', SMOOTH, ...
   'TRANSFORM',TRANSFORM,'EITRANSFORM',EITRANSFORM,'EGOAlg',EGOAlg,...
   'TolExpI',TolExpI,'pEst',pEst,'p0',p0,'pLow',pLow,'pUpp',pUpp,...
   'SAMPLEF',SAMPLEF,'KEPS',KEPS,'LCBB',LCBB,'GEIG',GEIG,...
   'globalSolver',globalSolver, 'localSolver',localSolver,...
   'LOCAL',LOCAL,...
   'GOMaxFunc',Prob.GO.MaxFunc,'GOMaxIter',Prob.GO.MaxIter,...
   'GOmaxFunc1',Prob.GO.maxFunc1,'GOmaxFunc2',Prob.GO.maxFunc2,...
   'GOmaxFunc3',Prob.GO.maxFunc3,...
   'backupSolver1',backupSolver1,'backupSolver2',backupSolver2);

% ----------------------
% cgolib initialization
% ----------------------
% Set default epsRank to epsRank
%cgolib(50, struct('DEFAULTEPSRANK', epsRank));
cgolib(50, struct('DEFAULTEPSRANK', epsRank,'PRINTLEVEL',1000,'PRINTERRORS',1000));

% Create cgolib objects
if isinf(nMax)
   cgolibnmax = MaxFunc;
else
   cgolibnmax = min(MaxFunc,nMax);
end

CGOLIB.setid = cgolib(100, d, cgolibnmax, 20);

% Create data column
cgolib(103, CGOLIB.setid, 1);

% Transformation table
% --------------------

% Create transformations
% List of transformations:
%  Column       Transformation   Source
%    0          Data
%    1          Median           0
%    2          Log              0
%    3          Log              1
%    4          Neg. log         0
%    5          Neg. log         1
%    6          Neg. reciprocal  0
%    7          Neg. reciprocal  1

cgolib(111, CGOLIB.setid, [1, 2, 2, 3, 3, 4, 4], [0, 0, 1, 0, 1, 0, 1]);
for i = 0:7
   CGOLIB.daceid(i+1) = cgolib(200, CGOLIB.setid, i);
   cgolib(250, CGOLIB.daceid(i+1), struct('SAMPLEF', SAMPLEF, 'GEIG', GEIG, ...
       'KEPS', KEPS, 'LCBB', LCBB));
end
% HKH NOTE BUG
% TRANSFORMATION 5 is not working properly for glb_prob #5 (REPLACE=1)

% Add the available data
cgolib(102, CGOLIB.setid, X, 0, F);



% -----------------------------------------
% Generate Prob structures for sub problems
% -----------------------------------------

% --------------------------------------------------------------
DACEProb = DACEGlobalProb(Prob, 'DACE', 'dace_f', [], [],...
   pEst, p0, pLow, pUpp, CGOLIB, globalSolver, localSolver, ...
   backupSolver1, backupSolver2, PriSub);
% --------------------------------------------------------------

%gProb = CGOGlobalProb(Prob, goName, glob_f, glob_g, glob_H,...
%   glob_c, glob_dc, glob_d2c, x_L, x_U, x_D, x_LL, x_UU, SCALE,...
%   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
%   backupSolver1, backupSolver2, PriSub)

%--------------------------------------------------------------

dc  = Prob.FUNCS.dc;
if isempty(dc)
   glob_dc = [];
else
   glob_dc = 'ego_dc';
end
EGOProb = CGOGlobalProb(Prob, 'EGO', 'ego_f', [], [],...
   'ego_c', glob_dc, [], x_L, x_U, x_D, x_LL, x_UU, SCALE,...
   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
   backupSolver1, backupSolver2, PriSub);

EGOProb.EGO.k          = d;
EGOProb.GO.ProbL.EGO.k = d;

% --------------------------------------------------
% Special input in EGOProb
% --------------------------------------------------
if DebugPriLev > 0
   EGOProb.optParam.IterPrint          = 1;
   EGOProb.PriLevOpt                   = DebugPriLev;
   EGOProb.GO.ProbL.optParam.IterPrint = 1;
   EGOProb.GO.ProbL.PriLevOpt          = DebugPriLev;
end
EGOProb.CGO.epsRank              = epsRank;
EGOProb.CGO.EITRANSFORM          = EITRANSFORM;
EGOProb.GO.ProbL.CGO.EITRANSFORM = EITRANSFORM;

% --------------------------------------------------
% ixDist (index for fMin) is used in constraints defined
% in rbf_c, rbf_dc, rbf_d2c.  Not currently used
% --------------------------------------------------
EGOProb.ixDist  = [];
% EGOProb.GO.ProbL.ixDist  = [];
DACEProb.ixDist = [];
% DACEProb.GO.ProbL.ixDist  = [];
% --------------------------------------------------

%--------------------------------------------------------------
% Define current best point
% (x_min, O_min, fMin) at index fIdx and flag Feasible
%--------------------------------------------------------------

z = Fpen;
% Set infeasible points as infinity before check on minium
z(Fpen-F >= 1E-14)=Inf;

% Set integer infeasible points as infinity before check on minium
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

% Check if transform value OK
if TRANSFORM > 1
   switch TRANSFORM
   case 2
    if any(F <= 0) 
       fprintf('TRANSFORM = %d not possible, min(F)=%f, use TRANSFORM = 0\n',...
                TRANSFORM,min(F));
       TRANSFORM = 0;
       Result.CGO.TRANSFORM = 0;
    end
   case 3
    if any(F >= 0) 
       fprintf('TRANSFORM = %d not possible, max(F)=%f, use TRANSFORM = 0\n',...
                TRANSFORM,max(F));
       TRANSFORM = 0;
       Result.CGO.TRANSFORM = 0;
    end
   case 4
    if ~(all(F > 0) | all(F <0)) | any(F==0)  
       fprintf('TRANSFORM = %d not OK, min(F)=%f, max(F)=%f, #F==0 %d.',...
                TRANSFORM,min(F),max(F),sum(F==0));
       fprintf(' Use TRANSFORM = 0\n')
       TRANSFORM = 0;
       Result.CGO.TRANSFORM = 0;
    end
   end
end


% This can maybe be removed later on:
DACEProb = ProbCheck(DACEProb,globalSolver,checkType('glb'),checkType('uc'));
%DACEProb.GO.ProbL = ProbCheck(DACEProb.GO.ProbL,DACEProb.GO.localSolver, ...
%checkType('con'),checkType('uc'));
if isempty(IntVars)
   EGOProb  = ProbCheck(EGOProb,globalSolver,checkType('glc'));
else
   EGOProb  = ProbCheck(EGOProb,globalSolver,checkType('minlp'));
end

% ------------------------------------
% Input parameters used in DACEFit
% ------------------------------------
DACEProb.CGO.globalSolver   = globalSolver;
DACEProb.CGO.localSolver    = localSolver;
DACEProb.PriLev             = 0;
% ------------------------------------
% For DACEFit: Make a Prob structure for estimation of one theta value
% ------------------------------------
DACEProb1          = DACEProb;
DACEProb1.N        = 1;
DACEProb1.pEst     = 0;
DACEProb1.p        = p0;
DACEProb1.x_L      = DACEProb.x_L(1);
DACEProb1.x_U      = DACEProb.x_U(1);
DACEProb1.x_0      = DACEProb1.x_L+1;
DACEProb1.x_min    = [];
DACEProb1.x_max    = [];
DACEProb1.GO.ProbL = [];
DACEProb1.thEst    = 1;
DACEProb.DACEProb1 = DACEProb1;
DACEProb.thEst     = d;

DACEProb.GO.ProbL.thEst = d;

% ------------------------------------
% Input parameters used in EGOSolve
% ------------------------------------
EGOProb.CGO.globalSolver    = globalSolver;
EGOProb.CGO.localSolver     = localSolver;
EGOProb.PriLev              = 0;

if PriLev > 1
   fprintf('globalSolver %s;',globalSolver);
   fprintf('localSolver %s;',localSolver);
   fprintf('GOlocalSolver %s;',GOlocalSolver);
   fprintf('backupSolver1 %s;',backupSolver1);
   fprintf('backupSolver2 %s.',backupSolver2);
   fprintf('\n');
end

% ------------------------------------
% Lipschitz initialization
% ------------------------------------

Dx     = tomsol(30,X,X);
LipLow =  inf;
LipUpp = -inf;
for i=1:n-1
   L = abs(F(i+1:n)-F(i))./Dx(i+1:n,i);
   LipLow = min(LipLow,min(L));
   LipUpp = max(LipUpp,max(L));
end
% Slower Alternative
%for i=1:n-1
%    [minDx, ixDx,Dx]  = cgolib(113, DACEProb.CGOLIB.setid, X(:,i));
%    L = abs(F(i+1:n)-F(i))./Dx(i+1:n);
%    LipLow = min(LipLow,min(L));
%    LipUpp = max(LipUpp,max(L));
%end

%--------------------------------------------------------------
% IterPrint at iteration 0
%--------------------------------------------------------------
Its  = [];
%xInf = NaN*ones(d,1);      %currently not used
%fInf = NaN;                %currently not used
if PriLev > 1 | IterPrint
   fprintf('Iter %3d n %3d ', 0, n);
   fprintf('nFunc %3d ', nFunc);
   tt = time([3 2 4 5]);
   if tt(4) < 10
      fprintf('%d/%d %d:0%1d ', tt);
   else
      fprintf('%d/%d %d:%2d ', tt);
   end
   fprintf('        ');   %fprintf('Step %2d ', modN);
   %fprintf('               ');
   if ~isempty(fGoal) & ~isinf(fGoal)
      fprintf(' fGoal %8.5f', fGoal);
   end
   if Feasible
      fprintf(' fMinF %11.8f', fMin);
   else
      fprintf(' fMinI %11.8f', fMin);
   end
   maxF = max(F);

   %fprintf(' max(F) %7.3g', maxF);
   %fprintf(' median(F) %7.3g', median(F));
   %fprintf(' range(F) %7.3g', maxF-fMin);
   %fprintf(' ||x_U-x_L|| %7.3g (pDist)', pDist);

   fprintf(' max(F) %9.5g', maxF);
   fprintf(' med(F) %9.5g', median(F));
   fprintf(' rng(F) %9.5g', maxF-fMin);

   if PriLev > 2
      fprintf('\n   ');
      fprintf('pDist %9.5g', pDist);
      fprintf(' LipU %f. ',LipUpp)
      fprintf('LipL %f. ',LipLow)
   else
      fprintf(' pDist %9.5g', pDist);
   end

   fprintf('\n   ');
   xprint(O_min,'xMin:',' %12.8f',8)
   if ~isempty(xOptS) & PriLev > 2
      xprint(Prob.x_opt,'   xOpt:',' %12.8f',8)
   end
end
%--------------------------------------------------------------

% -----------------------------------------------------------
%   cross-validation of the dace model to find best TRANSFORM
% -----------------------------------------------------------
if isempty(TRANSFORM)
   if PriLev > 1 | IterPrint
      fprintf('   Cross-validation of the DACE model\n');
   end

   [theta, p, fk] = DACEFit(DACEProb, IterPrint);

   [valid1,maxv1,sumv1] = crossval(DACEProb, PriLev);

   allpositive = all(F > 0);
   allnegative = all(F < 0);
   transfmap   = [0 0 4];

   if allpositive
      transfmap(2) = 2;
   elseif allnegative
      transfmap(2) = 3;
   else
      transfmap(2) = 0;
   end

   if(transfmap(2) ~= 0)
      DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(transfmap(2), REPLACE);
      DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
      [theta2, p2, fk2] = DACEFit(DACEProb, IterPrint);
      [valid2,maxv2,sumv2] = crossval(DACEProb, PriLev);
   else
      valid2 = 0; maxv2 = 0; sumv2 = 0;
   end

   % This transformation is only good if all values are less than 0 or all
   % values are greater than 0. No values can be 0.
   if allnegative | allpositive
      DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(transfmap(3), REPLACE);
      DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
      [theta3, p3, fk3] = DACEFit(DACEProb, IterPrint);
      [valid3,maxv3,sumv3] = crossval(DACEProb, PriLev);
   else
      valid3 = 0; maxv3 = inf; sumv3 = inf;
   end

   % This code picks a transformation using the following rules:
   %
   %  Only one of the valid transformations are picked.
   %  Among the valid transformationes, the one with the smallest sumv
   %  value is picked.

   validall = [valid1, valid2, valid3];
   sumall   = [sumv1, sumv2, sumv3];
   maxall   = [maxv1, maxv2, maxv3];

   valididx      = find(validall==1);

   % If there are no valid transformations. Make the transformation(s)
   % with lowest maxv valid.
   if(isempty(valididx))
      [tmp, valididx] = min(maxall);
   end

   [tmp, minidx] = min(sumall(valididx));
   minidx = minidx(1);

   ix = valididx(minidx);
   TRANSFORM = transfmap(ix);

   DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(TRANSFORM, REPLACE);
   DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;

   % Set best theta, p, fk
   if ix == 2
      theta = theta2;
      p     = p2;
      fk    = fk2;
   elseif ix == 3
      theta = theta3;
      p     = p3;
      fk    = fk3;
   end
   % Recompute
   %[theta, p, fk] = DACEFit(DACEProb, IterPrint);

   Result.CGO.TRANSFORM = TRANSFORM;
else
   DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(TRANSFORM, REPLACE);
   % Default transformation: No transformation.
   if(isempty(DACEProb.CGOLIB.TRANSFORM))
      DACEProb.CGOLIB.TRANSFORM = 0;
      TRANSFORM                 = 0;
      Result.CGO.TRANSFORM      = TRANSFORM;
   end
   DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
   theta = [];
end

% -------------------------------------------
% EGOProb input used in ego_f and DaceSurf_f:
% -------------------------------------------
EGOProb.EGO.CGOLIB          = DACEProb.CGOLIB;
EGOProb.GO.ProbL.EGO.CGOLIB = DACEProb.CGOLIB;

if PriLev > 1 | IterPrint
   fprintf('   ');
   fprintf('TRANSFORM %d. ', TRANSFORM);
   fprintf('EITRANSFORM %d. ', EITRANSFORM);
   fprintf('REPLACE %d ', REPLACE);
   if PriLev > 2
      fprintf('Internal TRANSFORM %d. ', DACEProb.CGOLIB.TRANSFORM);
   end
   fprintf('\n');
end


% *********************************************************
% *************** Initialize for MAIN ITERATION LOOP ***************
% *********************************************************

LowExp   = 0;

FLOWER    = fIdx(1);
fMinIter  = 0;

NOUPDATE  = 0;             % Used in test for Inform = 4
SAME1     = 0;             % Used in test for Inform = 5
SAME2     = 0;             % Used in test for Inform = 6
MaxCycle  = 5;             % Used in tests for Inform 4,5,6
O_pre     = inf*ones(d,1);

% -------------- Result saving -------------------
TESTDIG  = [];         % Struct with convergence info to known x_min

% -------------- Initialization -------------------
%ExitFlag = 1;   %overwritten before used
%ExitText = 'Maximal number of iterations reached'; %overwritten before used
TIME0    = Prob.TIME0; % Measure used CPU time
DIGIT    = 0.1;        % Convergence in "DIGIT" digits to known x_min
cpumax   = 0;          % Flag if max CPU reached
convflag = 0;          % Type of convergence

%Update   = 1;         %overwritten before used

% n     = Number of points used in the interpolation
% nFunc = Total number of costly function evaluations
% nInit = Number of points used for startup, now or previous run
% Iter  = Number of iterations to obtain new sample points
Iter     = 0;
% modN  = Type of search step in current iteration
modN     = 0;

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

% *********************************************************
% *************** START MAIN ITERATION LOOP ***************
% *********************************************************

while nFunc < MaxFunc & convflag == 0

   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   Iter = Iter + 1;

   % Set CGO field parameters in global structure EGOProb
   % HKH Not used now, maybe later for good initial tries
   % [X] = cgolib(109, DACEProb.CGOLIB.setid, n);
   % EGOProb.CGO.X         = X;
   % Get all points (is this really necessary? No, seems a waste)
   % % F = cgolib(110, DACEProb.CGOLIB.setid, n, 1, [], 0);
   % EGOProb.CGO.F         = F;
   % EGOProb.x_0           = x_min;

   % Set CGO field parameters and others in structure DACEProb
   DACEProb.Name          = ['DACE problem ' num2str(Iter)];
   DACEProb.GO.ProbL.Name = DACEProb.Name;


   % -----------------------------------------------------------------------
   %  Maximize Likelihood function to estimate (theta,p) for DACE surface
   % -----------------------------------------------------------------------

   if (Iter == 1 & isempty(theta)) | Iter > 1
      [theta, p, fk] = DACEFit(DACEProb, IterPrint);
   end
   % Interpolate DACE Surface with optimal theta and p.
   ID = DACEProb.CGOLIB.daceid(DACEProb.CGOLIB.TRANSFORM+1);
   cgolib(202, ID, theta, p);

   %NHQ plot
   if PLOTFLAG
      % Update Sampled Points in figure ProbFIG
      figure(ProbFIG)
      plot3(O(1,end),O(2,end),F(end),'k.','MarkerSize',15)


      % Update and plot current DACE-surface
      X = cgolib(109, DACEProb.CGOLIB.setid, n);
      DaceF = zeros(1,size(X,2));
      for k = 1:size(X,2)
         DaceF(k) = DaceSurf_f(X(:,k),EGOProb);
      end
      EGOProb.SCALE    = Prob.CGO.SCALE;
      EGOProb.FUNCS.f  = 'DaceSurf_f';
      EGOProb.probFile = Prob.probFile;
      clf(DaceFIG)
      plotProblem(EGOProb,X,DaceF,DaceFIG);
      pause(0.5)
      EGOProb.FUNCS.f = 'ego_f';
      fprintf(' modN = %3d\n',modN)
   end

   if modN >= 0

      % ----------------------------------------
      % Maximization of the expected improvement
      % ----------------------------------------
      EGOProb.Name          = ['EGO problem ' num2str(Iter)] ;
      EGOProb.GO.ProbL.Name = EGOProb.Name;
      %EGOProb.FUNCS.f = 'ego_f';
      EGOResult       = EGOSolve(EGOProb);
      xNew            = EGOResult.x_k(:,1);
      [ExpI, yMin]    = cgolib(204,DACEProb.CGOLIB.daceid(DACEProb.CGOLIB.TRANSFORM+1),xNew);
      % Optimal value of expected improvement (sign change from min to max)
      ExpI            = -ExpI;   %NHQ Är inte ExpI == -EGOResult.f_k  ???
      %OK, beror på EITRANSFORM. Sant om EITRANFROM = 0.
      %Hur/varför ändra denna??

      EGOResult.Prob.x_opt    = EGOResult.x_k';
      EGOResult.Prob.f_opt    = EGOResult.f_k;
      %EGOResult.Prob.probFile = Prob.probFile;
      %plotProblem(EGOResult.Prob)
      if PLOTFLAG
         clf(ExpFIG)
         plotProblem(EGOResult.Prob,X+0.0001,[],ExpFIG,1);
      end

   elseif modN == -1
      % -----------------------------------
      % Find global minimum of DACE surface
      % -----------------------------------
      EGOProb.Name             = ['DACE Surface minimum problem ' num2str(Iter)] ;
      EGOProb.FUNCS.f          = 'DaceSurf_f';
      EGOProb.GO.ProbL.Name    = EGOProb.Name;
      EGOProb.GO.ProbL.FUNCS.f = 'DaceSurf_f';
      EGOResult                = EGOSolve(EGOProb);
      xNew                     = EGOResult.x_k(:,1);

      if PLOTFLAG
         % Add DACE-surface minimum (feasible)
         figure(DaceFIG)
         plot3(EGOResult.x_k(1),EGOResult.x_k(2),EGOResult.f_k,'r*','MarkerSize',18)
      end

      EGOProb.FUNCS.f          = 'ego_f';
      EGOProb.GO.ProbL.FUNCS.f = 'ego_f';

   elseif modN == -2
      % ----------------------------------------
      % Maximization of the Kriging Variance
      % ----------------------------------------
      % NOTE - DOES NOT WORK BECAUSE mse_f must be rewritten to use cgolib
      % to compute
      EGOProb.Name             = ['Kriging Variance ' num2str(Iter)] ;
      EGOProb.FUNCS.f          = 'mse_f';
      EGOProb.GO.ProbL.Name    = EGOProb.Name;
      EGOProb.GO.ProbL.FUNCS.f = 'mse_f';
      EGOResult                = EGOSolve(EGOProb);
      xNew                     = EGOResult.x_k(:,1);
      EGOProb.FUNCS.f          = 'ego_f';
      EGOProb.GO.ProbL.FUNCS.f = 'ego_f';
   end

   if PLOTFLAG
      disp('done modN')
      %keyboard
   end

   % New point in original space
   if SCALE
      O_new   = tomsol(9, x_L, xNew, x_D);
   else
      O_new   = xNew;
   end


   % ---------------------------------------------------------
   % Clean TOMLAB global variables, because switch of problem
   % Compute surface errors in the new point xNew (O_new)
   % ---------------------------------------------------------
   NLP_x=[]; NLP_f=[]; NARG = [];

   EGOProb.FUNCS.f   = 'DaceSurf_f';
   DACENew           = nlp_f(xNew,EGOProb);
   if ~isempty(xOptS)
      snOptf         = nlp_f(xOptS,EGOProb);
      snOptg         = norm(nlp_g(xOptS,EGOProb));
      if isnan(snOptg)
         snOptH = [];
      else
         snOptH = nlp_H(xOptS,EGOProb);
      end
   end

   EGOProb.FUNCS.f = 'ego_f';
   NLP_x=[]; NLP_f=[]; NARG = [];

   % Distance between new point and best point found or closest point in X
   [onB, doX, doM] = statGN(EGOProb.x_L,EGOProb.x_U,xNew,x_min,X,epsX,IntVars,Reals);
   % Distance between new point and known optimal point, if given
   if ~isempty(xOptS)
      doO = min(tomsol(30,xNew,xOptS));
   else
      doO = [];
   end

   fRed = 0;

   % *************** UPDATE ***************
   % Remove information for response surface problem

   % Update  = -3; infeasible (xNew []), n > d+1, will remove worst infeasible point
   % Update  = -5; infeasible (xNew []), n = d+1, do nothing
   % Update  =  1; minDist > 1E-7. New point is distance 1E-7 from closest point
   % Update  =  2; 0 < minDist <= 1E-7, too close. Keep new, but remove closest point
   % Update  = -4; minDist == 0, n > d+1 (could remove a point, Update=-3,not now)
   % Update  = -4; minDist == 0, n = d+1

   if isempty(O_new) % Infeasibility problem
      fNew    = Inf;
      fPen    = Inf;
      SAME1   = 0;
      SAME2   = SAME2+1;
      % CGOLIBCHG
      [npoints] = cgolib(107, DACEProb.CGOLIB.setid);
      if npoints > d+1
         % It is possible to remove a point from X
         Update  = -3;
         minDist = -3;
      else
         Update  = -5;
         minDist = 0;
      end
   else
      if all(O_new == O_min), SAME1 = SAME1 + 1; else SAME1 = 0; end
      if all(O_new == O_pre), SAME2 = SAME2 + 1; else SAME2 = 0; end
      % CGOLIBCHG
      % Try to use functions from cgolib, not have multiple copies of X
      % Measure distance from xNew to X
      [minDist, ixDist,Dx] = cgolib(113, DACEProb.CGOLIB.setid, xNew);
      % ixDist is returned 0-based
      ixDist = ixDist + 1;

      if minDist ~= 0
         if minDist > 1E-7
            Update = 1;
         else
            Update = 2;
         end
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
         fRed  = fMin - fPen;
      else
         [npoints] = cgolib(107, DACEProb.CGOLIB.setid);
         if npoints > d+1
            % It is possible to remove a point from X
            % Update  = -2;
            % But should we do it, depends on the method in next step
            % Right now keep all points
            Update  = -4;
         else
            Update  = -4;
         end
         fNew    = NaN;
         fPen    = Inf;
      end
   end

   time  = fix(clock);

   if PriLev > 2
      fprintf('   Update %d ',Update);
      fprintf('minDist (to X set) %f ',minDist);
      fprintf('\n');
   elseif PriLev > 1 | IterPrint
      fprintf('\n');
   end

   if ~isnan(fNew)
      % Computed in call to cgolib(113,...)
      % Dx     = tomsol(30,X,xNew);
      L      = abs(F(1:n-1)-F(n))./Dx(1:n-1,1);
      LipLow = min(LipLow,min(L));
      LipUpp = max(LipUpp,max(L));
   end

   % HKH Avoid stop, must continue with new types of search steps in next Iter
   %if 0 & minDist == 0
   %   if PriLev > -1
   %      %xAlt
   %      fprintf('Convergence: New point identical to previous point');
   %      fprintf(' fMin %f,',fMin);
   %      fprintf(' %d iterations\n',Iter);
   %      fprintf('Maybe global optimum found?\n');
   %   end
   %   ExitFlag = 0;
   %   ExitText = 'Algorithm generates same points, maybe global optimum found';
   %   break;
   %end
   if Update == 1      % minDist > 1E-7
      % CGOLIBCHG
      F     = [F;fNew];
      Fpen  = [Fpen;fPen];
      F00   = [F00;f00New];
      if ~isempty(CcNew)
         Cc = [Cc,CcNew];
      end
      % X   = [X xNew]; % Saved in CGOLIB
      O     = [O O_new];
      O_pre = O_new;
      cgolib(102, DACEProb.CGOLIB.setid, xNew(:), 0, fNew);
      n     = n + 1;
      NOUPDATE = 0;
   elseif Update == 2  % 0 < minDist <= 1E-7
      if fPen < Fpen(ixDist)  % slight improvement
         % if Fpen(ixDist) > fNew always add new point, but delete close point
         % CGOLIBCHG
         % In cgolib this is done by adding the new point and then deleting the
         % old point, so that the new point will take the place of the old point.

         cgolib(102, DACEProb.CGOLIB.setid, xNew(:), 0, fNew);
         cgolib(104, DACEProb.CGOLIB.setid, 1, ixDist-1, 0);

         F(ixDist)       = fNew;
         Fpen(ixDist)    = fPen;
         F00(ixDist)     = f00New;
         if ~isempty(CcNew)
            Cc(:,ixDist) = CcNew;
         end
         %      X(:,ixDist) = xNew;
         O(:,ixDist)     = O_new;
         NOUPDATE        = 0;
      else
         % if Fpen(ixDist) < fNew skip adding new point.
         NOUPDATE = NOUPDATE+1;
      end
   elseif minDist == 0 % Update == -4 or -5
      NOUPDATE = NOUPDATE+1;
   else % Update == -3
      % Infeasible problem, ill-conditioning?
      % Remove most infeasible point
      [maxPen ixDist] = max(Fpen-F);
      if maxPen == 0
         % If all points feasible, remove point with largest f(x)
         [maxPen ixDist] = max(Fpen);
      end

      % CGOLIBCHG
      % Remove the point from cgolib. Use packing strategy 1.
      cgolib(104, DACEProb.CGOLIB.setid, 1, ixDist-1, 1);

      ix          =  [1:ixDist-1,ixDist+1:n];
      F           = F(ix);
      Fpen        = Fpen(ix);
      % X = X(:,ix);
      O           = O(:,ix);
      F00         = F00(ix);
      if ~isempty(CcNew)
         Cc       = Cc(:,ix);
      end
      n           = n - 1;
      if ixDist < fIdx
         fIdx     = fIdx -1;
      end
      %Update = 3;
      NOUPDATE = NOUPDATE+1;
   end

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

   % Is the function value in the new point less than fMin?
   % If feasible, compare fNew,fMin, otherwise fPen and fMin
   % Or now feasible, but previously not?
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
   end

   %if yMin >= 0
   %   yPred = (1-ExpI)*yMin;
   %else
   %   yPred = (1+ExpI)*yMin;
   %end

   yPred = yMin-ExpI;

   if any(TRANSFORM == [0 1])  % No transformation made.
      yNew = fNew;
      fPred = yPred;
   elseif TRANSFORM == 2  % log(y) transformation made.
      if fNew == 0
         % If the transformation is bad, change the transformation.
         DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(0, REPLACE);
         DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
         TRANSFORM = 0;
         yNew      = fNew;
         fPred     = yPred;
      else
         yNew      = log(fNew);
         fPred     = exp(yPred);
      end
   elseif TRANSFORM == 3  % -log(-y) transformation made.
      if fNew == 0
         % If the transformation is bad, change the transformation.
         DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(0, REPLACE);
         DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
         TRANSFORM = 0;
         yNew      = fNew;
         fPred     = yPred;
      else
         yNew       = -log(-fNew);
         fPred      = -exp(-yPred);
      end
   elseif TRANSFORM == 4  % -1/y transformation made.
      if fNew == 0
         % If the transformation is bad, change the transformation.
         DACEProb.CGOLIB.TRANSFORM          = GetCGOLIBTransform(0, REPLACE);
         DACEProb.GO.ProbL.CGOLIB.TRANSFORM = DACEProb.CGOLIB.TRANSFORM;
         TRANSFORM = 0;
         yNew      = fNew;
         fPred     = yPred;
      else
         yNew      = -1/fNew;
         fPred     = -1/yPred;
      end
   end
   % Compute relative surface error as
   % (Actual value - Predicted surface value) / Actual value
   surfErr = yNew-DACENew;
   if yNew ~= 0, surfErr = surfErr/abs(yNew); end

   %  yMin=cgolib(110,DACEProb.CGOLIB.setid,1,1,fIdx,DACEProb.CGOLIB.TRANSFORM);

   yMin = cgolib(209, DACEProb.CGOLIB.daceid(DACEProb.CGOLIB.TRANSFORM+1));

   if isempty(fGoal) | isinf(fGoal)
      RelErr = NaN;
   else
      if fGoal == 0
         RelErr = abs(fMin-fGoal);
      else
         RelErr = abs(fMin-fGoal)/abs(fGoal);
      end
   end

   % Save information from the current iteration
   Its.Iter(Iter)             = Iter;
   Its.n(Iter)                = n;
   Its.nFunc(Iter)            = nFunc;
   Its.modN(Iter)             = modN;
   Its.fMin(Iter)             = fMin;
   Its.fRed(Iter)             = fRed;
   Its.RelErr(Iter)           = RelErr;
   Its.surfErr(Iter)          = surfErr;
   Its.FLOWER(Iter)           = FLOWER;
   Its.fMinIter(Iter)         = fMinIter;
   Its.fNew(Iter)             = fNew;
   Its.onBxNew(Iter)          = onB;
   Its.distxNew2X(Iter,:)     = doX;
   Its.distxNew2xMin(Iter,:)  = doM;
   if 1 & ~isempty(xOptS)
      if isnan(snOptg)
         snOptE = NaN;
      else
         snOptE = sum(eig(snOptH)<0);
      end
      % dXO              = min(tomsol(30,xOptS,X));
      dXO                = cgolib(113, DACEProb.CGOLIB.setid, xOptS);

      Its.snOptf(Iter)   = snOptf;
      Its.snOptg(Iter)   = snOptg;
      Its.snOptE(Iter)   = snOptE;
      Its.dXO(Iter)      = dXO;
   end

   % Lipschitz variables
   Its.LipUpp(Iter)      = LipUpp;
   Its.LipLow(Iter)      = LipLow;

   if SAVE == 1
      fMinIdx  = fIdx(1);
      rngState = rand('state');
      % Get the rows
      [X] = cgolib(109, DACEProb.CGOLIB.setid, n);
      saveCGO(1,'ego',Name,O,F,X,F_m,F00,Cc,nInit,Fpen,fMinIdx,rngState);
   end

   if PriLev > 1 | IterPrint
      fprintf('Iter %3d n %3d ', Iter, n);
      fprintf('nFunc %3d ', nFunc);
      fprintf('%d/%d %d:%d ', time([3 2 4 5]));
      fprintf('Step %2d ', modN);
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
      % fprintf('at %3d fNew %11.8f ', FLOWER, fNew);
      if ~isnan(RelErr)
         fprintf(' RelE%8.5g ', RelErr);
      end
      fprintf(' DErr %10.6g', surfErr);
      %NOTfprintf(' fLoc %11.8f', min_sn);

      fprintf('\n');
      % Special EGO printout
      fprintf('   ');
      if ~isempty(CcMin)
         fprintf('CcMin:  ');
         fprintf('%11.8f ', CcMin);
      end
      if fk < -1E200
         fprintf('fE %8.5f ', -inf);
      else
         fprintf('fE %8.5f ',-fk);
      end
      if modN >= 0
         fprintf('fPred %8.5f ', fPred);
         fprintf('I%%');
         fprintf(' %7.4g ',100*ExpI);
         fprintf('Q%%');
         if yMin ~= 0, fprintf(' %7.4g ',100*abs(ExpI/yMin)); end
      end
      if PriLev > 2
         fprintf('[%d] ',onB);
         if isempty(IntVars)
            fprintf('doX %f ',doX);
            fprintf('doM %f ',doM);
         else
            fprintf('doX %f/%d ',doX);
            fprintf('doM %f/%d ',doM);
         end
         if ~isempty(xOptS)
            fprintf('doO %f ',doO);
         end
         fprintf('fRed %f. ',fRed)
         fprintf('LipU %f. ',LipUpp)
         fprintf('LipL %f. ',LipLow)
      end
      if TRANSFORM ~= 0
         if modN >= 0
            fprintf('yMin %11.8f yNew %11.8f ', yMin, yNew);
            fprintf('yPred %8.5f ', yPred);
         else
            fprintf('yNew %11.8f',yNew);
         end
      end
      fprintf('\n');
      xprint(O_new,'   xNew:',' %12.8f',10)
      if d > 5
         xprint(    p,'      p:',' %10.6f',12)
         xprint(theta,'  theta:',' %10.6f',12)
      else
         xprint([theta,p],'   th,p:',' %8.4f',10);
      end
      if ~isempty(xOptS) & PriLev > 2
         fprintf('   dXO %f ',dXO);
         fprintf('DACE(xOpt) %8.5f ',snOptf);
         fprintf('||DACE-g(xOpt)|| %f ',snOptg);
         fprintf('Negeig(DACE-H(xOpt)) %d ',snOptE);
      end

   end

   % ********** CONVERGENCE TEST **********

   if n >= nMax, convflag = 7; end
   if convflag == 0
      convflag = CGOisClose(fGoal,fMin,fTol,nFunc,Iter,PriLev);
   end
   % HKH Is any of these convergence tests redundant?
   % Good value for MaxCycle?
   if convflag == 0
      if NOUPDATE > MaxCycle
         convflag = 4;
      elseif SAME1 > MaxCycle
         convflag = 5;
      elseif SAME2 > MaxCycle
         convflag = 6;
      end
   end
   % Convergence test on scaled expected improvement
   if modN >= 0
      if fMin == 0
         ExpII = ExpI;
      else
         ExpII = ExpI/abs(fMin);
      end
      if ExpII <= TolExpI & fRed <= 0
         LowExp = LowExp + 1;
         if LowExp >= 3
            convflag = 10;
            ExitText = 'Expected improvement low for three iterations';
            if PriLev > 1 | IterPrint
               fprintf('\n');
               fprintf('   Convergence: ');
               fprintf('ExpI %e ',ExpI);
               fprintf('ExpI/|fMin| %e <',ExpII);
               fprintf(' TolExpI %e',TolExpI);
               fprintf('. Iteration %d',Iter);
               fprintf('\n');
               fprintf('   %s\n',ExitText);
            end
            break;
         else
            if PriLev > 1 | IterPrint
               fprintf('\n');
               fprintf('   LowExp increased %d: ',LowExp);
               fprintf('ExpI %e ',ExpI);
               fprintf('ExpI/|fMin| %e <',ExpII);
               fprintf(' TolExpI %e',TolExpI);
               fprintf('. Iteration %d',Iter);
               fprintf('\n');
            end
         end
      elseif ExpII <= TolExpI & fRed > 0
         if PriLev > 1 | IterPrint
            fprintf('ExpI %e ',ExpI);
            fprintf('ExpI/|fMin| %e <',ExpII);
            fprintf(' TolExpI %e',TolExpI);
            fprintf(', BUT fRed %f >0\n',fRed);
         end
         yRed = fRed;
         if fMin~=0, yRed = yRed/abs(fMin); end
         if yRed > TolExpI
            % Low ExpI, but improvement is too high, reset LowExp
            LowExp = 0;
         end
      elseif LowExp > 0
         yRed = fRed;
         if fMin~=0, yRed = yRed/abs(fMin); end
         if yRed > TolExpI
            % Too high improvement and ExpI value, reset LowExp
            LowExp = 0;
         end
      end
   elseif fRed > 0
      yRed = fRed;
      if fMin~=0, yRed = yRed/abs(fMin); end
      if yRed > TolExpI
         % Too high improvement in other type of step, reset LowExp
         LowExp = 0;
      end
   end
   if PriLev > 0 & LowExp > 0
      fprintf('----- >>> LowExp %d\n',LowExp);
   end

   % ----------------------------------
   % Initial p value for next iteration
   % ----------------------------------
   if pEst == 0
      DACEProb.x_0 = log(theta);
   elseif pEst == 1
      DACEProb.x_0 = [log(theta);p(1)];
   else
      DACEProb.x_0 = [log(theta);p];
   end
   % disp( [DACEProb.x_L, DACEProb.x_0, DACEProb.x_U])

   % ------------------------------------------
   % Print digits of convergence to given fGoal
   % ------------------------------------------
   while RelErr <= DIGIT & DIGIT > 1E-6
      k = abs(log10(DIGIT));
      if (IterPrint | PriLev > 0) %& convflag == 0
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


   % ------------------------------------------
   % Determine what type of step next iteration
   % ------------------------------------------
   if EGOAlg == 1
      % Use EXPI when progress, switch to 1 step DACE if no progress
      if modN >= 0 & fRed <= 0
         modN = -1;
      else
         modN = modN + 1;
      end
   elseif EGOAlg == 2
      % Always use EXPI, except after LowExp == 2
      if modN >= 0 & fRed <= 0
         if LowExp > 1
            modN = -1;
         else
            modN = LowExp;
         end
      else
         %modN = modN + 1;
         modN = 0;
      end
   end

end % Main iteration loop

% *******************************************************
% *************** END MAIN ITERATION LOOP ***************
% *******************************************************


% Get the rows
[X] = cgolib(109, DACEProb.CGOLIB.setid, n);

% X transformed to original coordinates already saved in O matrix.


%NHQ
if RMSError & d <= 4 
   % Calculate RMS error of the interpolated surface.
   % Performed for d = 2,3 or 4.
   % Crossvalidate surface before setting free cgolib.
   fprintf(' Starting valuation of the final RBF-surface.\n')
   % Always use original values in order to compare with real f(x).
   DACEProb.CGOLIB.TRANSFORM = 0;
   DACEProb.GO.ProbL.CGOLIB.TRANSFORM = 0;
   % Fit surface using ML and Interpolate DACE Surface.
   [theta,p,fk] = DACEFit(DACEProb, IterPrint); % theta,p,fk never used???
   % Interpolate DACE Surface with optimal theta and p.
   ID = DACEProb.CGOLIB.daceid(DACEProb.CGOLIB.TRANSFORM+1);
   cgolib(202, ID, theta, p);
   NLP_x=[]; NLP_f=[]; NARG = [];

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
            sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - cgolib(205,0,xD);
         end
      end
      if PLOTFLAG
         figure('Name','Difference EGO Surface vs. True Function');
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
               sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - cgolib(205,0,xD);
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
                  sFx(kk) = eval([Prob.FUNCS.f '(xF,Prob)']) - cgolib(205,0,xD);
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


% Free the cgolib set and all the surfaces connected to it.
cgolib(101, DACEProb.CGOLIB.setid);

% SAVE RESULTS
fMinIdx  = fIdx(1);
rngState = rand('state');
% Create struct output with warmstart information; save results on cgoSave.mat
Result.CGO.WarmStartInfo = saveCGO(2,'ego', Name,O,F,X,F_m,F00,Cc,...
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
   fprintf('constraints are violated')
   if isempty(Result.x_k)
      % HKH FIX - something is wrong
      disp(SCALE)
      disp(fMin)
      disp(Fpen)
      disp(Fpen-fMin)
      xprinti(find(Fpen-fMin >  0),'>0:')
      xprinti(find(Fpen-fMin == 0),'=0:')
      xprinti(find(Fpen-fMin <  0),'<0:')
      [fxxx,ixxx] = min(Fpen);
      x_k  = X(:,ixxx);
      if SCALE
         Result.x_k      = tomsol(9,x_L,x_k,x_D);
      else
         Result.x_k      = x_k;
      end
   end
end

% if SCALE
%   Result.x_k      = tomsol(9,x_L,X(:,find(F==fMin)),x_D);
% else
%   Result.x_k      = X(:,find(F==fMin));
% end
% if isempty(Result.x_k)
%   if SCALE
%     Result.x_k      = tomsol(9,x_L,X(:,find(Fpen==fMin)),x_D);
%   else
%     Result.x_k      = X(:,find(Fpen==fMin));
%   end
%   fprintf('Warning: Optimal point is found to be infeasible: ');
%   fprintf('constraints are violated')
% end

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
%cMin = zeros(max([size(Prob.c_L) size(Prob.c_U)]),size(Result.x_k,2));
%
%if dCon > 0
%   for i=1:size(Result.x_k,2)
%      %cMin = [cMin,nlp_c(Result.x_k(:,i), Prob, varargin{:})];
%      cMin(:,i) = nlp_c(Result.x_k(:,i), Prob, varargin{:});
%      nCon = nCon + 1;
%   end
%end
%if dLin > 0 & size(Result.x_k,2)==1
%   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at best x_k
%end
if dLin > 0
   Result.Ax    = Prob.A*Result.x_k;  % Linear constraint value at all x_k
end
Result.CGO.DACEProb = DACEProb;
Result.CGO.EGOProb  = EGOProb;
Result.CGO.fGoal    = fGoal;
Result.CGO.Its      = Its;
Result.Iter         = Iter;     % Number of iterations
Result.FuncEv       = nFunc;
Result.GradEv       = -1;
Result.HessEv       = -1;
Result.ConJacEv     = -1;
Result.ExitFlag     = 0;
Result.ExitText     = ['Tried ' num2str(nFunc) ' f(x), using ' ...
   num2str(n) ', startup ' num2str(nInit), ' ',ExDText];
Result.SolverAlgorithm = [Result.SolverAlgorithm ...
   '. Global solver ' globalSolver  ...
   '. Local solver ' localSolver];
if cpumax & convflag == 0
   Result.Inform   = 9;
   Result.ExitText = [Result.ExitText '. Max CPU reached. '];
elseif convflag == 7
   Result.Inform   = convflag;
   Result.ExitText = [Result.ExitText, '. All feasible integers tried'];
else
   Result.Inform   = convflag;
end
Result.DIGIT      = TESTDIG;

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

Result            = endSolve(Prob,Result);

% -------------- End of main EGO routine -------------------
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
function EGOResult = EGOSolve(EGOProb)
%-----------------------------------------------------------
if nargin < 1
   error('EGOSolve needs one parameter, the structure EGOProb');
end

% Pick up local parameters
globalSolver = EGOProb.CGO.globalSolver;
%localSolver  = EGOProb.CGO.localSolver;    %currently not used
PriLev       = EGOProb.PriLev;
EGOProb.xInit = max(100,4*EGOProb.N);

%[EGOProb.x_L,EGOProb.x_0,EGOProb.x_U]
EGOResult           = tomRun(globalSolver,EGOProb,PriLev-4);
%EGOResult           = tomRun(globalSolver,EGOProb,2);
%EGOResult           = tomRun('multiMin',EGOProb,2);

%ExitFlag = EGOResult.ExitFlag;

%-----------------------------------------------------------
function [theta, p, fk] = DACEFit(DACEProb, IterPrint)
%-----------------------------------------------------------

% Pick up local parameters
globalSolver = DACEProb.CGO.globalSolver;
localSolver  = DACEProb.CGO.localSolver;
PriLev       = DACEProb.PriLev;
d            = DACEProb.d;

% Compute the transformation to be used
cgolib(112, DACEProb.CGOLIB.setid, DACEProb.CGOLIB.TRANSFORM);

%[DACEProb.x_L,DACEProb.x_0,DACEProb.x_U]

% x_0 0 or 1 will give immediate stop

% Solve 1-dim problem in theta
DACEProb1    = DACEProb.DACEProb1;
DACEResult1  = tomRun(localSolver, DACEProb1, PriLev-4);
% DACEResult1  = tomRun(localSolver, DACEProb1, 2);

pEst = DACEProb.pEst;
x_L  = DACEProb.x_L;
x_U  = DACEProb.x_U;
x_D  = x_U-x_L;
if pEst > 0
   DACEProb.pEst = 0;
   x_0           = DACEProb.x_0;
   DACEProb.N    = d;
   DACEProb.x_L  = x_L(1:d);
   DACEProb.x_U  = x_U(1:d);
   DACEProb.p    = DACEProb1.p;
end

DACEProb.x_0       = DACEResult1.x_k*ones(d,1);
% Adjust lower bounds if 1-dim solution give value less than x_L(2:d)
if any(DACEResult1.x_k < x_L(1:d))
   Adjust             = 1;
   % xprint(DACEProb.x_L,'xL0')
   x_L(2:d)           = DACEProb.x_L(2:d)-2;
   DACEProb.x_L(2:d)  = x_L(2:d);
   % xprint(DACEProb.x_L,'xL1')
else
   Adjust             = 0;
end

% Solve full problem in theta, with x_0 = 1-dim solution
DACEResult2 = tomRun(localSolver, DACEProb, PriLev-4);
%DACEResult2 = tomRun(localSolver, DACEProb, 2);
% xprint(DACEResult2.Prob.x_0,'x_0')
x_k = DACEResult2.x_k;

if Adjust
   % Reset lower bounds
   x_L(2:d)           = min(x_L(2:d),x_k(2:d));
end

if pEst > 0
   DACEProb.pEst      = pEst;
   DACEProb.x_0       = x_0;
   DACEProb.x_0(1:d)  = x_k;
   % DACEProb.x_0(d+1)  = DACEProb.p;
   % DACEProb.x_0(d+1)  = 1.9;
   DACEProb.N         = length(x_U);
   DACEProb.x_L       = x_L;
   DACEProb.x_U       = x_U;
   DACEResult2        = tomRun(localSolver, DACEProb, PriLev-4);
   %DACEResult2        = tomRun(localSolver, DACEProb, 2);
else
   DACEProb.x_L(2:d)  = x_L(2:d);
   DACEProb.x_0(1:d)  = x_k;
end

%DACEProb.xInit = min(100,20*DACEProb.N);
%DACEResult4 = tomRun(globalSolver, DACEProb, 2);

% Use multiMin for a limited amount of searches
% But limit the upper bounds to avoid initial values at extreme crap points
M              = max(20,DACEProb.N);
if pEst > 0
   DACEProb.xInit = daceInit(M,ceil(M/2),[x_L(1:d);0.9*x_0(d+1:end)],[x_U(1:d)-3;x_U(d+1:end)],[],0.05*norm(x_D));
else
   DACEProb.xInit = daceInit(M,ceil(M/2),x_L(1:d),x_U(1:d)-3, [],0.05*norm(x_D));
end
% DACEProb.PriLevOpt = 5;
DACEResult3    = tomRun(globalSolver, DACEProb, PriLev-4);
%DACEResult3    = tomRun(globalSolver, DACEProb, 2);
% xprint(DACEResult3.multiMin.fOpt,'fOp3')

%M              = max(20,2*DACEProb.N);
%DACEProb.xInit = daceInit(M,ceil(M/2),x_L(1:d),x_U(1:d)-3, [],0.05*norm(x_D));
%DACEProb.xInit = -max(20,2*DACEProb.N);
%DACEResult3 = tomRun(globalSolver, DACEProb, PriLev-4);
%DACEResult4    = tomRun(globalSolver, DACEProb, 2);
%xprint(DACEResult4.multiMin.fOpt,'fOp4')


fk3      = DACEResult3.f_k;
fk2      = DACEResult2.f_k;
fk1      = DACEResult1.f_k;

if IterPrint > 1
   fprintf('   Theta solutions: ');
   fprintf('fk1 %13.7f. ',fk1);
   fprintf('fk2 %13.7f. ',fk2);
   fprintf('fk3 %13.7f. ',fk3);
   fprintf('FuncEv:');
   fprintf(' %4d',DACEResult1.FuncEv);
   fprintf(' %4d',DACEResult2.FuncEv);
   fprintf(' %6d. ',DACEResult3.FuncEv);
   fprintf('Inform: ');
   fprintf(' %d ',DACEResult1.Inform);
   fprintf(' %d ',DACEResult2.Inform);
   fprintf(' %d ',DACEResult3.Inform);
   fprintf('\n');
end

if fk1 <= fk2 & fk1 <= fk3
   x_k      = DACEResult1.x_k*ones(d,1);
   fk       = fk1;
   if pEst > 0
      x_k   = [x_k;DACEProb1.p];
   end
elseif fk2 <= fk3
   x_k      = DACEResult2.x_k(:,1);
   fk       = fk2;
else
   x_k      = DACEResult3.x_k(:,1);
   fk       = fk3;
end
theta    = exp(x_k(1:d));

%if DACEResult2.f_k < DACEResult3.f_k
%	keyboard
%end
%ExitFlag = DACEResult.ExitFlag;    %currently not used

if length(x_k) > d
   p = x_k(d+1:end);
   if isscalar(p)
      p = p*ones(d,1);
   end
   %if PriLev >= 0
   %   xprint(p,'   pNew:  ');
   %end
else
   p = DACEProb.p;
end

%-----------------------------------------------------------
% function [valid, maxv, sumv] = crossval(DACEProb, PriLev)
%-----------------------------------------------------------

function [valid, maxv, sumv] = crossval(DACEProb, PriLev)

[valid, maxv, sumv] = cgolib(206, DACEProb.CGOLIB.daceid(DACEProb.CGOLIB.TRANSFORM+1));

if PriLev > 1
   fprintf('Internal TRANSFORM %d ',DACEProb.CGOLIB.TRANSFORM);
   fprintf('maxv %d ',maxv);
   fprintf('sumv %d ',sumv);
   fprintf('valid %d',valid)
   fprintf('\n');
end


% % --------------------------------------------------------------
% function [invR, detR, pRank] = getinvR(R, epsRank);
% % --------------------------------------------------------------
% n          = size(R,1);
% [U S V]    = svd(R);
% S_inv      = zeros(n,1);
% S11        = S(1,1);
% detR       = S11;
% S_inv(1) = 1/S11;
% for i = 2:n
%   Sii = S(i,i);
%   if Sii > epsRank*S11
%     pRank    = i;
%     detR     = detR*Sii;
%     S_inv(i) = 1/Sii;
%   else
%     break;
%   end
% end
% invR = V(:,1:pRank) * diag(S_inv(1:pRank)) * U(:,1:pRank)';

% fprintf('invR: pRank %d n %d detR %20.4e\n',pRank,n,detR);

% --------------------------------------------------------------
function [CGOTRANSFORM] = GetCGOLIBTransform(TRANSFORM, REPLACE)
% --------------------------------------------------------------
%
% TRANSFORM must be passed non-empty and within the ranges,
% otherwise GetCGOLIBTransform will return empty.
%
% The computations here is based on the "Transformation table" above.

CGOTRANSFORM = [];

if(isempty(REPLACE))
   REPLACE = 0;
end

if(~isempty(TRANSFORM))
   % The transformation has to be within the range.
   if(TRANSFORM < 0 | TRANSFORM > 4)
      return
      %TRANSFORM = [];
   else
      % If the user explicitly has chosen the median transformation,
      % set REPLACE to 1 instead, and let TRANSFORM be 0.
      if(TRANSFORM == 1)
         TRANSFORM = 0;
         REPLACE   = 1;
      end
      % Now TRANSFORM can take one of the values: {0 2 3 4}. Not 1,
      % the median. If median transformation is used, then REPLACE ==
      % 1. We want TRANSFORM to take one of the values {0 1 2 3}
      % instead, in order to map it to the transformation columns.
      if(TRANSFORM ~= 0)
         TRANSFORM = TRANSFORM-1;
      end
      % Finally, check the REPLACE value. Should we apply the median
      % transformation?
      if(REPLACE == 1)
         CGOTRANSFORM = TRANSFORM*2+1;
      else
         CGOTRANSFORM = TRANSFORM*2;
      end
   end
end

% --------------------------------------------------------------
function [TRANSFORM, EITRANSFORM, EGOAlg, TolExpI, pEst, p0, pLow, pUpp, ...
        SAMPLEF, KEPS, LCBB, GEIG] = ...
    getEGOProbVars(Prob)
%getEGOProbVars(Prob, x_L, x_U)
% --------------------------------------------------------------

if isempty(Prob.CGO)
   TRANSFORM    = []; EITRANSFORM  = []; EGOAlg      = []; TolExpI  = [];
   pEst         = []; p0           = []; pLow        = []; pUpp     = [];
   SAMPLEF      = []; KEPS         = []; LCBB        = []; GEIG     = [];
else
   if isfield(Prob.CGO,'TRANSFORM')
      TRANSFORM = Prob.CGO.TRANSFORM;
   else
      TRANSFORM = [];
   end
   if isfield(Prob.CGO,'EITRANSFORM')
      EITRANSFORM = Prob.CGO.EITRANSFORM;
   else
      EITRANSFORM = [];
   end
   if isfield(Prob.CGO,'EGOAlg')
      EGOAlg = Prob.CGO.EGOAlg;
   else
      EGOAlg = [];
   end
   if isfield(Prob.CGO,'TolExpI')
      TolExpI = Prob.CGO.TolExpI;
   else
      TolExpI = [];
   end
   if isfield(Prob.CGO,'pEst')
      pEst = Prob.CGO.pEst;
   else
      pEst = [];
   end
   if isfield(Prob.CGO,'p0')
      p0 = Prob.CGO.p0;
   else
      p0 = [];
   end
   if isfield(Prob.CGO,'pLow')
      pLow = Prob.CGO.pLow;
   else
      pLow = [];
   end
   if isfield(Prob.CGO,'pUpp')
      pUpp = Prob.CGO.pUpp;
   else
      pUpp = [];
   end
   if isfield(Prob.CGO,'SAMPLEF')
       SAMPLEF = Prob.CGO.SAMPLEF;
   else
       SAMPLEF = [];
   end
   if isfield(Prob.CGO,'KEPS')
       KEPS = Prob.CGO.KEPS;
   else
       KEPS = [];
   end
   if isfield(Prob.CGO,'LCBB')
       LCBB = Prob.CGO.LCBB;
   else
       LCBB = [];
   end
   if isfield(Prob.CGO,'GEIG')
       GEIG = Prob.CGO.GEIG;
   else
       GEIG = [];
   end
end


%NOT, keep empty,if isempty(TRANSFORM),    TRANSFORM = -1; end
if isempty(EITRANSFORM),  EITRANSFORM = 1; end
EITRANSFORM = min(2, max(0, EITRANSFORM));

if isempty(EGOAlg),       EGOAlg = 1; end
if isempty(TolExpI),      TolExpI = 1E-7; end

% Set default p0, pLow and pUpp dependent on pEst
if isempty(pEst),         pEst = 0; end

pEst = max(0,pEst);

if pEst == 0
   if isempty(p0),   p0   = 1.99; end
elseif pEst == 1
   if isempty(p0),   p0   = 1.99; end
   if isempty(pLow), pLow = 1; end
   if isempty(pUpp), pUpp = 1.9999; end
else
   d = Prob.N;
   if isempty(p0),   p0   = 1.99*ones(d,1); end
   if isempty(pLow), pLow = ones(d,1); end
   if isempty(pUpp), pUpp = 1.9999*ones(d,1); end
end

%************************************
% The rest is just a saving of code
%************************************

% % contour plot of the Likelihood function
%
% u = linspace(1E-5,0.05,50);
% v = linspace(1E-5,0.005,50);
% for j=1:length(u)
%    j
%    for i=1:length(v)
%       x=[u(j);v(i)];
%       ff(i,j)=DACE_f(x,DACEProb);
%       %ff(i,j)=feval(DACEProb.FUNCS.f,x,DACEProb);
%   end
%end
%contour(u,v,ff,50);
%save 'vars' u v ff
%end
%
% % contour plot of model function
%
%o = ones(n-1,1);
%y = y(1:length(y)-1);
%u = linspace(x_L(1),x_U(1),100);
%v = linspace(x_L(2),x_U(2),100);
% %for j=1:length(u)
%   for i=1:length(v)
%      x=[u(j);v(i)];
%      r = zeros(n-1,1);
%      for kk = 1:n-1
%         r(kk) = exp(-( sum( theta.*(abs( x-X(:,kk) )).^p) ));
%      end
%      %ff(i,j)=nlp_f(x,Prob, varargin{:});
%      ff(i,j) = my_hat + r'*invR*(y-my_hat*o);
%   end
%end
%contour(u,v,ff,75);
%end
%
%plot(DACEResult.GLOBAL.C(1,:),DACEResult.GLOBAL.C(2,:),'*');


% MODIFICATION LOG
%
% 981127  mbk  DACEProb is defined without dace_prob.m.
%              DACEProb.x_0 is set in glb_prob.
% 981127  hkh  Add calls to ini/endSolve, globalSave/Get.
% 981208  hkh  Init of yMin = Inf; xBest= Prob.x_0;
%              Add Result.ExitFlag. Fix error return of Result.
% 981208  mbk  Init of yMin = Inf; xBest= Prob.x_0; ABORTED
%              Return the smallest function value of those initially
%              computed and corresponding x-value.
% 001118  hkh  Major revision
% 010304  hkh  Use tomSolve
% 010718  hkh  Use glbFast, revised similar to glbSolve
% 010726  hkh  Change local / solver handling and printing.
% 010903  hkh  General initialization of optPar vector to -999.
% 020103  hkh  Major revision taking many elements from rbfSolve
% 020616  hkh  Use input TolExpI and TRANSFORM.
% 020616  hkh  Fix bug calling SOL solvers, estimating internal derivatives
% 020820  hkh  Print sub problem results dependent on PriLev
% 020820  hkh  Try local search for all possible global solutions
% 020823  hkh  Error in transforming linear constraints if SCALE = 1
% 020824  hkh  Accept local solution only if linear constraints fulfilled
% 020828  hkh  Use GetSolver to select default local solver
% 021020  hkh  Check dCon > 0 before using xnargin
% 021020  hkh  Use Fpen=F+sum(max(0,dev)) to handle infeasible initial points
% 030222  hkh  Add check for empty Result.x_k, use best infeasible solution
% 030514  hkh  Avoid size as variable name in randomtest, use p instead
% 030514  hkh  Use rbf_c interface if pEst set, because x is expanded
% 030514  hkh  Also do xnargin on rbf_d2c (normally not used)
% 040104  hkh  Use NARGSAVE to set NARG, must set nlp_c=nlp_xc = [];
% 040104  hkh  sum2 changed to sumv2, sum3 to sumv3
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin in DACEProb and EGOProb
% 040309  hkh  Revised functions gutmann and corners, added new sampleInts
% 040309  hkh  Set field Prob.MIP in DACEProb.MIP and EGOProb.MIP
% 040309  hkh  Set globalSolver='glcCluster' if Prob.MIP.IntVars is non-empty
% 040309  hkh  Making EGO handle mixed-integer problems
% 040309  hkh  Add CGO.RandState, to initialize the random generator
% 040309  hkh  Check for cycle of same points
% 040309  hkh  Set Result.Inform to convflag
% 040309  hkh  Revise feasibility handling, print fMinI, fMinF, use Feasible
% 040309  hkh  Special treatment of pure IP, with glcSolve and glcFast
% 040310  hkh  Add new CGO.pEst parameter, if true estimate || ||_p (default)
% 040323  hkh  Revised with similar initial point features as rbfSolve
% 040323  hkh  Set cTol and bTol in suboptimizations, avoiding infeasibility
% 040323  hkh  Default for Percent if dLin+dCon > 0, use -5000 instead of -d
% 040404  hkh  Default for pEst=0, for SCALE = 1
% 040404  hkh  Major revision, global solvers do not return empty
% 040406  hkh  Save fMinIdx (=fIdx(1)) on cgoSave.mat
% 040414  hkh  Do not set MENU or LargeScale in EGOProb and DACEProb
% 040414  hkh  Return Result.c_k and Result.Ax
% 040415  hkh  Safe guard if localSolver return > 1 solutions
% 040415  hkh  Avoid LOCAL=1 if pure IP problems
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 040923  hkh  Wrong XL XU when pEst=0, DACEProb wrong.
% 040923  hkh  Changed text to "Cannot proceed", not Failure, if ill-conditioned
% 040929  hkh  sampleInts, x_L(j), offset from 0, must be added
% 041005  hkh  Add localSolver and globalSolver to SolverAlgorithm output
% 041123  hkh  Change call to tomRun
% 050218  hkh  New function daceInts, round experimental design points
% 050218  hkh  Add check on duplicate points in sampleInts and daceInts
% 050218  hkh  Use daceInts instead of sampleInts when for DACE strategies
% 050218  hkh  Add new output field Result.CGO with tuning params used
% 050218  hkh  Made new getInvR using SVD to compute safe inv(R) (too slow)
% 050321  hkh  Must define n  = size(X,2); for integer case in initial strategy
% 050322  hkh  Wrong check for duplicates in daceInts, also round only IntVars
% 050322  hkh  Pure IP: Test if all integer values tried (nMax == n)
% 050322  hkh  Test if all integer values tried, if all continuous fixed
% 050322  hkh  Add test if all variables are fixed
% 050427  hkh  Change call to daceInit for new version
% 050504  hkh  TolExp=1E-9, instead of 1E-6. Fixed bugs calling daceInit
% 050504  hkh  Use Prob.optParam.eps_Rank for epsRank, 1E-11 instead of 1E-14
% 050504  hkh  Test epsRank   = max(1E-14,epsRank);
% 060104  frhe DACEProb no longer a copy of main Prob.
% 060106  frhe Transformation on ego_f added. Controlled with EITRANSFORM.
% 060106  frhe DACEProb.dDim is equal to num vars of problem.
% 060207  hkh  Changed glcFast to glcDirect as default, added glcDirect setup
% 060207  hkh  Add new input nSample
% 060207  hkh  Use new function expDesign for initial experimental design
% 060608  hkh  Print TRANSFORM, EITRANSFORM and REPLACE if IterPrint | PriLev>0
% 060608  hkh  Add EITRANSFORM as Result.CGO output
% 060630  hkh  Make a structure of Result.Dig3 with fields Iter,FuncEv,CPUtime
% 060630  hkh  Add DIRECT methods as initial design methods, Percent > 100
% 060814  med  FUNCS used for callbacks instead
% 060816  hkh  Prob.MIP.IntVars also set if isempty(Prob.MIP), expDesign crash
% 060817  hkh  Test if finite bounds. Define Result before tests
% 060818  hkh  Avoid picking up Prob.f_Low, not used now
% 061004  ango Safe writing of cgoSave.mat. Added Result.CGO.WarmStartInfo
% 061124  hkh  Lower letter in case tests
% 070221  hkh  Add output field Result.DIGIT for 1,...,6 digits convergence
% 070222  hkh  Revise IntVars handling, use new format
% 070531  hkh  DACEProb simple bounds wrongly defined as matrix if pEst = 1
% 070531  hkh  Change minDist tolerance check from 1E-5 to 1E-7
% 071006  hkh  Do init of random generator only if no warm start, NaN no init
% 071007  hkh  Save random generator state in rngState, reinit rand if WarmStart
% 080115  hkh  Major revision, make similar with with arbfMIP and rbfSolve
% 080225  hkh  End Major revision
% 080414  hkh  Prob.GO input in Result.CGO, change CGOGLobalProb call
% 080414  hkh  Add nTrial, CLHMethod to Result.CGO
% 080417  hkh  Use relative constraint violation in fPen computation
% 080617  frhe Made REPLACE>1 transformation monotonic
% 080617  frhe log10 used in REPLACE>1 in accordance with help
% 080630  hkh  Use cgo_fc instead of nlp_f for costly f(x) (and costly c(x))
% 080711  hkh  Use routine saveCGO for warm start info saving for all solvers
% 080711  hkh  Save warm start info every iteration in cgoSave.mat (cgoSave1)
% 080711  hkh  Avoid estimation of Hessian in endSolve setting HessEv = -1
% 081104  hkh  Avoid crash in eig(H) when DaceSurf_g is bad with NaN values
% 081105  hkh  Avoid extra lines of output
% 081105  hkh  Make getProbVars solver dependent
% 081105  hkh  Add input/output nTrial in call to getProbVars
% 090425  hkh  Print theta,p instead of p,theta, and less decimals
% 090425  hkh  Print fE = -inf if fE <-1E300
% 090425  hkh  Set pUpp default as 1.9999, not 2
% 090902  hkh  if TRANSFORM>0, check if F is OK, otherwise set TRANSFORM=0
% 090909  hkh  cgolibnmax=min(MaxFunc,nMax); avoid cgolib crash for big nMax
% 090919  hkh  Avoid solving same DACEFit problem 3 times, now only once
% 090919  hkh  Do 1 instead of 5 DACE fitting solutions at start, cgolib(202)
% 091001  hkh  Add IterPrint to DACEfit input, print theta information
% 100222  hkh  Major revision of ExpI convergence handling, also change TolExpI=1E-7.
% 100222  hkh  Print exit information only if PriLev >1 or IterPrint
% 100222  hkh  Only use new close point (||< 1E-7) if fPen < F(ixDist), otherwise no update
% 100222  hkh  Revised DACEFit to handle pEst > 0
