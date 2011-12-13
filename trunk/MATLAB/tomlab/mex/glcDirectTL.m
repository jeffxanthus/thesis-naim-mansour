% glcDirect implements a modified version of the algorithm DIRECT.
% 
% glcDirect solves general constrained mixed integer global optimization
%
% glcDirect implements the DIRECT algorithm by Donald R. Jones 
% presented in the paper "DIRECT", Encyclopedia of Optimization, 
% Kluwer Academic Publishers, 2001.
% The algorithm is expanded to handle nonlinear and linear equalities, and
% linear inequalities.
%
% glcDirect solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Recommendation: Put the integers as the first variables !
% Put low range integers before large range integers
% Linear constraints are specially treated
% Equality constraints are added as penalties to the objective.
% Weights are computed automatically, assuimg f(x) scaled to be roughly 1
% at optimum. Otherwise scale f(x)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"
%
% Calling syntax:
%
% function Result = glcDirectTL(Prob)
%
% INPUT PARAMETERS
%
% Prob    Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say GLCF
%   FUNCS.c   The routine to compute the nonlinear constraint, say GLCC
%             A call to tomFiles.m or glcAssign.m sets these fields. 
%   x_L       Lower bounds for each element in x. Try to make a tight bound
%   x_U       Upper bounds for each element in x. Try to make a tight bound
%
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%
%   PriLevOpt Print Level. This controls both regular printing from
%             glcDirect and the amount of iteration log information
%             to print. 
%
%              0 = Silent
%              1 = Warnings and errors printed. Iteration log
%                  on iterations improving function value.
%              2 = Iteration log on all iterations.
%              3 = Log for each function evaluation.
%              4 = Print list of parameter settings.
%
%             See optParam.IterPrint for more information on
%             iteration log printing.
%
%   WarmStart If true, >0, glcDirect reads the output from the last
%             run from Prob.glcDirect.WarmStartInfo if it
%             exists. If it doesn't exist, glcDirect attempts to
%             open and read warm start data from mat-file
%             glcDirectSave.mat. glcDirect uses this warm start
%             information to continue from the last run.
%
%   MaxCPU    Maximum CPU Time (in seconds) to be used. 
%
% glcDirect   Structure with DIRECT algorithm specific parameters. Fields used:
%
%   fcALL     =0 (Default). If linear constraints cannot be feasible anywhere
%             inside rectangle, skip f(x) and c(x) computation for middle point
%             =1. Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Do not update rates of change 
%             for the constraints.
%             =2. Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Update rates of change
%             constraints.
%
%   useRoC    = 1 (Default). Use original Rate of Change (RoC) for constraints
%             to weight the constraint violations in selecting which rectangle
%             divide.
%             = 0. Avoid RoC, giving equal weights to all constraint violations.
%             Suggested if difficulty to find feasible points. For problems
%             where linear constraints have been added among the nonlinear
%             (NOT RECOMMENDED; AVOID!!!), then option useRoc=0 has been 
%             successful, whereas useRoC competely fails. 
%             = 2. Avoid RoC for linear constraints, giving weight one to these
%             constraint violations, whereas the nonlinear constraints use RoC.
%             = 3. Use RoC for nonlinear constraints, but linear constraints are
%             not used to determine which rectangle to use.
%
%   BRANCH    =0. Divide rectangle by selecting the longest side, if ties use
%             the lowest index. This is the Jones DIRECT paper strategy.
%             =1. First branch the integer variables, selecting the variable
%             with the least splits. If all integer variables are split, split
%             on the continuous variables as in BRANCH=0.  DEFAULT!
%             Normally much more efficient than =0 for mixed-integer problems
%             =2. First branch the integer variables with 1,2 or 3 possible values,
%             e.g [0,1],[0,2] variables, selecting the variable with least splits. 
%             Then branch the other integer variables, selecting the variable
%             with the least splits. If all integer variables are split, split
%             on the continuous variables as in BRANCH=0. 
%             =3. Like =2, but use priorities on the variables, similar to
%             mipSolve, see Prob.MIP.VarWeight.
%
%   RECTIE    When minimizing the measure to find which new rectangle to try to 
%             get feasible, there are often ties, several rectangles have the same
%             minimum. RECTIE = 0 or 1 seems reasonable choices. Rectangles with low
%             index are often larger then the rectangles with higher index.
%             Selecting one of each type could help, but often =0 is fastest.
%             = 0 Use the rectangle with value a, with lowest index (original)
%             = 1 (Default): Use 1 of the smallest and 1 of largest rectangles
%             = 2 Use the last rectangle with the same value a, not the 1st
%             = 3 Use one of the smallest rectangles with same value a
%             = 4 Use all rectangles with the same value a, not just the 1st
%
%   EqConFac  Weight factor for equality constraints when adding to objective
%             function f(x) (Default value 10). The weight is computed as
%             EqConFac/"right or left hand side constant value", e.g. if
%             the constraint is Ax <= b, the weight is EqConFac/b
%             If DIRECT just is pushing down the f(x) value instead of
%             fulfilling the equality constraints, increase EqConFac
%
%   AxFeas    Set nonzero to make glcDirect skip f(x) evaluations, when the
%             linear constraints are infeasible, and still no feasible point
%             has been found.  The default is 0. Value 1 demands fcALL == 0.
%             This option could save some time if f(x) is a bit costly, however
%             overall performance could on some problems be dramatically worse
%           
%   fEqual    All points with function values within tolerance fEqual are
%             considered to be global minima and returned. Default 1E-10
%
%   LinWeight RateOfChange = LinWeight*|a(i,:)| for linear constraints.
%             Balance between linear and nonlinear constraints.  Default 0.1. 
%             The higher value, the less influence from linear constraints
%
%   alpha     Exponential forgetting factor in RoC computation, default 0.9
%
%   AvIter    How many values to use in startup of RoC computation before
%             switching to exponential smoothing with forgetting factor alpha
%             Default 50.
%
% optParam    Structure in Prob, Prob.optParam. 
%             Defines optimization parameters. Fields used:
%   MaxFunc   Maximal number of function evaluations, default 10000
%   MaxIter   Maximal number of iterations, default 200
%   IterPrint Print iteration #, # of evaluated points, best f(x) and x for 
%             each iteration, where the best point has improved
%   bTol      Linear constraint feasibility tolerance
%   cTol      Nonlinear constraint feasibility tolerance
%   fGoal     Goal for function value, if empty not used
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%   eps_x     Convergence tolerance in x. All possible rectangles are 
%             less than this tolerance (scaled to (0,1) )
%             See the output field maxTri.
%   EpsGlob   Global/local weight parameter, default 1E-4
%
% MIP         Structure in Prob, Prob.MIP.
%             Defines integer optimization parameters. Fields used:
%   IntVars:  
%               If empty, all variables are assumed non-integer 
%               If islogical(IntVars) (=all elements are 0/1), then
%               1 = integer variable, 0 = continuous variable.
%               If any element >1, IntVars is the indices for integer variables
%
%             It is advised to number the integer values as the first
%             variables, before the continuous. The tree search will then
%             be done more efficiently.
%
%   fIP       An upper bound on the optimal f(x) value.  If empty, set as Inf.
%   xIP       The x-values giving the fIP value. 
%             If fIP empty and xIP given, fIP will be computed
%             if xIP nonempty, its feasibility is checked
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  c_k      Nonlinear constraints values at x_k
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ConstrEv Number of constraint evaluations(=FuncEv if nonlinear constraints)
%  maxTri   Maximum size of any triangle
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag 0 = Normal termination, max number of iterations /func.evals reached
%           2 = Some upper bounds below lower bounds
%           4 = Numerical trouble, and cannot continue
%           7 = Reached maxFunc or maxIter, NOT feasible
%           8 = Empty domain for integer variables
%           10= Input errors
%  Inform   1 = Function value f is less than fGoal
%           2 = Absolute function value f is less than fTol, only if fGoal = 0
%            or Relative error in function value f is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           3 = Maximum number of iterations done
%           4 = Maximum number of function evaluations done
%           5 = Maximum number of function evaluations would most likely
%               be too many in the next iteration, save warm start info, stop
%           6 = Maximum number of function evaluations would most likely
%               be too many in the next iteration, because
%               2*sLen >= maxFDim - nFunc, save warm start info, stop
%           7 = Space is dense
%           8 = Either space is dense, or MIP is dense
%           10= No progress in this run, return solution from previous one
%           91= Infeasible
%           92= No rectangles to work on
%           93= sLen = 0, no feasible integer solution exists
%           94= All variables are fixed
%           95= There exist free constraints
%
%  glcDirect  Substructure for glcDirect specific result data.
%   WarmStartInfo
%             Structure with warm start information. Use
%             WarmDefGLOBAL to reuse this information in another
%             run.
%   convFlag  Converge status flag from solver.
%
% To make a warm start possible, glcDirect saves the following
% information in the structure Result.glcDirect.WarmStartInfo and
% file glcDirectSave.mat (for internal solver use only):
%   C         Matrix with all rectangle centerpoints, in [0,1]-space.
%   D         Vector with distances from centerpoint to the vertices.
%   F         Vector with function values.
%   G         Matrix with constraint values for each point.
%   nIter     Number of iterations
%   Name      Name of the problem. Used for security if doing warm start
%   Split     Split(i,j) = # splits along dimension i of rectangle j
%   Tr        Tr(i) is the number of times rectangle i has been trisected.
%   fMinIdx   Indices of the currently best points.
%   fMinEQ    sum(abs(infeasibilities)) for minimum points, 0 if no equalities
%   glcfMin   Best function value found at a feasible point.
%   feasible  Flag indicating if a feasible point has been found.
%   ign       Rectangles to be ignored in the rect. selection procedure.
%   roc       Rate of change s, for each constraint
%   s0        Sum of observed rate of change s0 in the objective
%   t         t(i) is the total # splits along dimension i.
%
% USAGE:
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
%    Prob   = glcAssign('GLCF',x_L,x_U,'GLCF Test',A,b_L,b_U,'GLCC',c_L,c_U);
%    Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%    % Default values are now set for PriLevOpt, and structure optParam
%    % To change a value, examples are shown on the two next lines
%    Prob.optParam.MaxFunc = 500; % Change max number of function evaluations 
%    Prob.optParam.MaxIter = 200; % Change the number of iterations to 200
%
%    If there are integer variables, they may be set as additional input
%    to glcAssign, or directly as the next line shows:
%    Prob.MIP.IntVars = [1 3];  % 1st and third variables are integers
%
% -- Calling glcDirect
%
% Direct solver call:
%    Result = glcDirectTL(Prob);
%    PrintResult(Result);
%           
% Driver call, including printing with level 2:
%      Result = tomRun('glcDirect',Prob,2);
%           
% -- Defining an objective function and nonlinear constraints
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
%
% -- Warm start
%
% To make a restart (warm start), just set the warm start flag, and call
% glcDirect once again:
%
%      Prob.WarmStart = 1;
%      Result         = tomRun('glcDirect', Prob, 2);
%
% glcDirect will read warm start information from the glcDirectSave.mat file
%
% Another warm start (with same MaxFunc) is made by just calling tomRun again
%      Result = tomRun('glcDirect', Prob, 2);
%
% To make a restart from the warm start information in the Result
% structure, make a call to WarmDefGLOBAL before calling glcDirect.
% WarmDefGLOBAL moves information from the Result structure
% to the Prob structure and sets the warm start flag, Prob.WarmStart = 1; 
%
%      Prob = WarmDefGLOBAL('glcDirect', Prob, Result);
% where Result is the result structure returned by the previous run.
% A warm start (with same MaxIter) is done by just calling tomRun again
%
%      Result = tomRun('glcDirect', Prob, 2);
%
% To make another warm start (with same MaxFunc), repeat the two lines:
%      Prob   = WarmDefGLOBAL('glcDirect', Prob, Result);
%      Result = tomRun('glcDirect', Prob, 2);
%
%
% IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM:
%
% The DIRECT algorithm only reaches the variable bounds in the limit.
% Therefore convergence for global optimum where components are on the bounds
% is slow. 
% One remedy is to reduce lower bounds with a tolerance, say 1E-4, and add
% a similar tolerance 1E-4 to the upper bounds that might be reached.
% Another possibility is to fix a variable on its bound by setting the lower
% and upper bounds equal.
%
% The bound problem could be even worse for constrained problems. If a 
% nonlinear problem only is feasible when a certain component is on its bound,
% then DIRECT will never get feasible in reasonable time. This problem could
% be even worse if the problem and the constraint is mixed-integer.
% The solution is to fix variables on their bounds, or, at least temporary,
% set the troublesome variables as integers (if the bounds are integer-valued)
% As for unconstrained problems, relaxing the bounds slightly, might also be
% a way to go.
%
% Also try to avoid linear equality constraints, especially if some of the
% variables are integers. If say, two linear constraints are equalities, it
% is always possible to eliminate two variables, and reduce the dimension of
% the problem. 
%
% Always try to reduce the dimension as much as possible when using the
% DIRECT algorithm, and try to shrink the box, defined by the lower and
% Put low range integers before large range integers
% Linear constraints are specially treated
% Equality constraints are added as penalties to the objective.
% Weights are computed automatically, assuimg f(x) scaled to be roughly 1
% at optimum. Otherwise scale f(x)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Feb 15, 1999.   Last modified Oct 2, 2009.

function Result = glcDirectTL(Prob)

if nargin < 1
   error('glcDirect needs input structure Prob');
end

solvType=checkType('glc');

Prob=ProbCheck(Prob,'glcDirect',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = DefPar(Prob,'MaxCPU',2000000000);

% Pick up input parameters from the Prob structure:
x_L = Prob.x_L(:);   % Lower bounds
x_U = Prob.x_U(:);   % Upper bounds

if isempty(x_L) | isempty(x_U) | any(isinf(x_L) | isinf(x_U))
   Prob = preSolve(Prob);
   x_L = Prob.x_L;
   x_U = Prob.x_U;
end

n  = max(length(x_L),length(x_U));   % Problem dimension

% Check for Inf and set to lower values.
x_L(isinf(x_L)) = -10000;
x_U(isinf(x_U)) =  10000;

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
      error('glcDirect: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
% glcDirect mex takes IntVars as vector of indices
IntVars = find(IV);

isMIP     = ~isempty(IntVars);

PriLev    = Prob.PriLevOpt;          % Print level
MaxFunc   = Prob.optParam.MaxFunc;   % Maximum number of function evaluations
MaxIter   = Prob.optParam.MaxIter;   % Number of iterations
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).

fTol      = Prob.optParam.eps_f;     % Relative tolerance for fGoal
EpsGlob   = Prob.optParam.EpsGlob;   % global/local weight parameter. 
bTol      = Prob.optParam.bTol;      % Linear constraints feasibility tolerance
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance
xTol      = Prob.optParam.eps_x;     % Tolerance for rectangle sizes
                                     % (scaled to (0,1) )

% Safeguard
if MaxIter < 0
   MaxIter = [];
end
if MaxFunc < 0
   MaxFunc = [];
end

if isinf(fGoal) | isempty(fGoal),   fGoal = -1E300; end

Result                 = ResultDef(Prob);
Result.Solver          = 'glcDirect';
Result.SolverAlgorithm = 'Constrained DIRECT - C Implementation';


%A   = Prob.A;        % Linear constraint matrix
b_L = Prob.b_L(:);   % Lower bounds, linear constraints
b_U = Prob.b_U(:);   % Upper bounds, linear constraints
c_L = Prob.c_L(:);   % Lower bounds, nonlinear constraints
c_U = Prob.c_U(:);   % Upper bounds, nonlinear constraints

mNL = max(length(c_L), length(c_U));
c_L = [c_L(:) zeros(1, mNL-length(c_L))];
c_U = [c_U(:) zeros(1, mNL-length(c_U))];

mL = max(length(b_L), length(b_U));
b_L = [b_L(:) zeros(1, mL-length(b_L))];
b_U = [b_U(:) zeros(1, mL-length(b_U))];

fcn  = Prob.FUNCS.f;
cfcn = Prob.FUNCS.c;
if isempty(cfcn)
   if xnargin(fcn) > 1 
      probPtr = Prob;
   else
      probPtr = [];
   end
else
   n1 = xnargin(fcn);
   n2 = xnargin(cfcn);
   if (n1 == 2 & n2 == 2) | ( n1 == 3 & n2 == 3) | (n1 == 3 & n2 == 2)
      % Assume varargin if n1 == 3,  n2 == 2
      probPtr = Prob;
   elseif n1 == 1 & n2 == 1
      probPtr = [];
   else
      if(isa(fcn,'function_handle'))
         fprintf('Objective  function %s has %d input arguments\n',func2str(fcn),n1);
      else
         fprintf('Objective  function %s has %d input arguments\n',fcn,n1);
      end
      if(isa(cfcn,'function_handle'))
         fprintf('Constraint function %s has %d input arguments\n',func2str(cfcn),n2);
      else
         fprintf('Constraint function %s has %d input arguments\n',cfcn,n2);
      end
      error('glcDirect: Same number of input arguments is needed');
   end
end

if isfield(Prob.MIP,'fIP')
   fIP = Prob.MIP.fIP;
else
   fIP = []; % No estimate of optimal solution
end
if isfield(Prob.MIP,'xIP')
   xIP = Prob.MIP.xIP;
else
   xIP = []; % No solution estimate of optimal solution
end
if isempty(fIP)
   if ~isempty(xIP)
      % User has supplied an x
      if length(xIP) ~= n
         fprintf('Input Prob.MIP.xIP has wrong length %d\n',length(xIP));
         error('Illegal input Prob.MIP.xIP')
      end
      if isMIP
         xIP(IntVars) = round(xIP(IntVars));
      end
      fIP = nlp_f(xIP, Prob);  % Function value at x
   else
      fIP = Inf;
   end
elseif isMIP
   if ~isempty(xIP)
      xIP(IntVars) = round(xIP(IntVars));
   end
end

% Check if xIP is feasible
if ~isinf(fIP) & ~isempty(xIP)
   if ~isempty(Prob.A)
      Ax = Prob.A*xIP;
      if any(xIP < x_L | xIP > x_U) | ...
         any(b_L-Ax > bTol) | any(Ax-b_U > bTol)
         if IterPrint | PriLev > 0
            disp('Input Prob.MIP.xIP is not linear feasible');
            disp('Reject value of Prob.MIP.fIP');
         end
         fIP = Inf;     
         xIP = [];     
      end
   end
end

% glcDirect parameters
GD = DefPar(Prob,'glcDirect',[]);

options = DefPar(GD, 'options', []);

if Prob.WarmStart == 1
   % Restart with values from previous run.

   Name1 = deblank(Prob.Name);  % Problem name
   
   WarmStartInfo = DefPar(GD, 'WarmStartInfo', []);
   
   % If the WarmStart structure exists, then use that warm start data.
   if ~isempty(WarmStartInfo)
     Name = WarmStartInfo.Name;
     if strcmp(Name1,Name)
       if(PriLev >= 2)
         fprintf(['Restart with %d sampled points from last ' ...
                  'runs.\nUsing warm start data from structure: ' ...
                  'Prob.glcDirect.WarmStartInfo.\n'],...
                  WarmStartInfo.nFunc);
       end
     else
       Prob.WarmStart = 0;
       if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf(['Impossible to do restart using data from ' ...
                  'structure: Prob.glcDirect.WarmStartInfo.\n']);
       end
     end

   % If no WarmStart structure exists, look for a file.
   else
     filename = 'glcDirectSave.mat';
     
     if(exist(filename,'file')~=2)
       fprintf(['Couldn''t find warm start info structure nor warm ' ...
                'start file in path.\nNo warm start possible.\n']);
       Prob.WarmStart = 0;
     else
       load(filename,'Name');

       if strcmp(Name1,Name)
         WarmStartInfo = load(filename,'nFunc','nIter','feasible', ...
                              'fMinIdx','s0','glcfMin','fMinEQ', ...
                              'C', 'F','D','G','roc','Tr','ign', ...
                              'iL','iU','Split','t','nFuncTot',...
                              'nRectTot');
        
         if PriLev >= 2
           fprintf(['Restart with %d sampled points from last ' ...
                    'runs.\nUsing warm start data from file: %s.\n'], ...
                   WarmStartInfo.nFunc, filename);
         end
         
       else
         Prob.WarmStart = 0;
         if PriLev >= -1000
           fprintf('Previous run was with Problem %s\n',Name);
           fprintf('This run is with Problem %s\n',Name1);
           fprintf(['Impossible to do restart using data from file: ' ...
                    '%s.\n'], filename);
         end
       end
     end
   end
end

if Prob.WarmStart == 1 
  nIter0    = WarmStartInfo.nIter;
  nFuncTot0 = WarmStartInfo.nFuncTot;
else
  nIter0        = 0;
  nFuncTot0     = 0;
  WarmStartInfo = [];
end

AxFeas    = DefPar(GD,'AxFeas',0);
fEqual    = DefPar(GD,'fEqual',1e-10);
fcALL     = DefPar(GD,'fcALL',0);
useRoC    = DefPar(GD,'useRoC',1);
BRANCH    = DefPar(GD,'BRANCH',1);
RECTIE    = DefPar(GD,'RECTIE',1);
EqConFac  = DefPar(GD,'EqConFac',10);
LinWeight = DefPar(GD,'LinWeight',0.1);
alpha     = DefPar(GD,'alpha',0.9);
AvIter    = DefPar(GD,'AvIter',50);

AxFeas    = min(fcALL,AxFeas);  % AxFeas must be 0 if fcALL is 0

if ~isempty(AxFeas)
  options.AXFEAS = DefPar(options, 'AXFEAS', AxFeas);
end

if ~isempty(PriLev)
  options.PRILEV = DefPar(options, 'PRILEV', PriLev);
end

if ~isempty(MaxFunc)
  options.MAXFUNC = DefPar(options, 'MAXFUNC', MaxFunc);
end

if ~isempty(MaxIter)
  options.MAXITER = DefPar(options, 'MAXITER', MaxIter);
end

if ~isempty(AvIter)
  options.AVITER = DefPar(options, 'AVITER', AvIter);
end

if ~isempty(Prob.WarmStart)
  options.WARMSTART = DefPar(options, 'WARMSTART', double(Prob.WarmStart > 0));
end

if ~isempty(IterPrint)
  options.ITERPRINT = DefPar(options, 'ITERPRINT', IterPrint);
end

if ~isempty(MaxCPU)
  options.MAXCPU = DefPar(options, 'MAXCPU', MaxCPU);
  options.MAXCPU = max(0, min(options.MAXCPU, 2000000000));
end

if ~isempty(fTol)
  options.FUNTOL = DefPar(options, 'FUNTOL', fTol);
end

if ~isempty(xTol)
  options.VARTOL = DefPar(options, 'VARTOL', xTol);
end

if ~isempty(EpsGlob)
  options.GLWEIGHT = DefPar(options, 'GLWEIGHT', EpsGlob);
end

if ~isempty(cTol)
  options.NLCONTOL = DefPar(options, 'NLCONTOL', cTol);
end

if ~isempty(bTol)
  options.LCONTOL = DefPar(options, 'LCONTOL', bTol);
end

if ~isempty(alpha)
  options.ALPHA = DefPar(options, 'ALPHA', alpha);
end

if ~isempty(LinWeight)
  options.LINWEIGHT = DefPar(options, 'LINWEIGHT', LinWeight);
end

if ~isempty(fEqual)
  options.FEQUAL = DefPar(options, 'FEQUAL', fEqual);
end

if ~isempty(EqConFac)
  options.EQCONFAC = DefPar(options, 'EQCONFAC', EqConFac);
end

if ~isempty(fcALL)
  options.FCALL = DefPar(options, 'FCALL', fcALL);
end

if ~isempty(useRoC)
  options.USEROC = DefPar(options, 'USEROC', useRoC);
end

if ~isempty(BRANCH)
  options.BRANCH = DefPar(options, 'BRANCH', BRANCH);
end

if ~isempty(RECTIE)
  options.RECTIE = DefPar(options, 'RECTIE', RECTIE);
end

if ~isempty(fIP)
  options.FIP = DefPar(options, 'FIP', fIP);
end

varWeight = DefPar(Prob.MIP, 'VarWeight');
Name = deblank(Prob.Name);  % Problem name

[convFlag, xMin, glcfMin, nFuncTot, nRectTot, nIter, minSplit, ...
 nFunc, WarmStartInfoOut, minPnts] = ...
    glcDirect(fcn, x_L, x_U, IntVars, varWeight, ...
              full(Prob.A), b_L, b_U, cfcn, ...
              Prob.c_L, Prob.c_U, ...
              fGoal, options, WarmStartInfo, probPtr);

% Add name to WarmStartInfo structure. It is TL-specific.
WarmStartInfoOut.Name = Name;
Result.glcDirect.WarmStartInfo = WarmStartInfoOut;
clear WarmStartInfoOut;

switch convFlag
 case 3
   ExitFlag = 0;
   Inform   = 3;
   if Prob.WarmStart
      Result.ExitText = ['Max iterations. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(nIter) ' iter.'];
   else
      Result.ExitText = 'Maximal number of iterations reached.';
   end
case 9
   ExitFlag = 0;
   Inform   = 9;
   if Prob.WarmStart
      Result.ExitText = ['CPU time limit reached. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(nIter) ' iter.'];
   else
      Result.ExitText = 'CPU time limit reached.';
   end
 
   % 4 and 5 are removed.
 case {2,4,5}
   ExitFlag = 0;
   if convFlag == 4
      Inform   = 5;
   elseif convFlag == 5
      Inform   = 6;
   else
      Inform   = 4;
   end
   if Prob.WarmStart
      Result.ExitText = ['Max f(x) evals. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(nIter) ' iter.'];
   else
      Result.ExitText = 'Maximum number of f(x) evaluations reached.';
   end
 case 1
   ExitFlag = 0;
   Inform   = 1;
   if Prob.WarmStart
      Result.ExitText = ['Converged to f(x) < fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(nIter) ' iter.'];
   else
      Result.ExitText = 'Converged to f(x) < fGoal.';
   end
 case 11
   ExitFlag = 0;
   Inform   = 2;
   if Prob.WarmStart
      Result.ExitText = ['Converged to fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(nIter) ' iter.'];
   else
      Result.ExitText = 'Converged to fGoal.';
   end
 case -1
   ExitFlag = 4;
   Inform   = 94;
   Result.ExitText = 'All variables are fixed.';
 case -2
   ExitFlag = 0;
   Inform   = 92;
   Result.ExitText = 'No rectangles left in search area.';
 case -4 % I guess this will never occur. -2 will occur instead.
   ExitFlag = 8;
   Inform   = 93;
   Result.ExitText = 'No feasible integer solution exists.';
 case {-5, -6}
   ExitFlag = 4;
   Inform   = 95;
   Result.ExitText = 'There exist free constraints.';
 case -7
  % Should never occur here. Test is done before call to glcDirect.
   ExitFlag = 10;
   Inform   = -10000;
   Result.ExitText = 'Invalid format on IntVars array.';
 otherwise
   ExitFlag = 4;
   Inform   = 99;
   Result.ExitText = 'Numerical trouble - EXIT!';
end

if Inform < 90
  if minSplit == -1
    Inform   = 7;
    Result.ExitText = 'Space is dense.';
  end
  if minSplit == -2
    % ExitFlag = 0; Does this only occur if convFlag = -1
    Inform   = 8;
    Result.ExitText = 'Either space is dense, or MIP is dense.';
  end
end

Result.ExitText = sprintf(['%s\nTried %i rectangles, computing %i ' ...
                          'f(x) values.'], Result.ExitText, nRectTot, ...
                          nFuncTot);

if minSplit < 0
   Result.maxTri   = 0;            % Maximal size of possible triangles
else
   Result.maxTri   = 3.^-minSplit; % Maximal size of possible triangles
end
Result.Iter     = nIter-nIter0;   % Number of iterations this run
Result.FuncEv   = nFuncTot-nFuncTot0;

Result.glcDirect.convFlag = convFlag;

feasible = Result.glcDirect.WarmStartInfo.feasible;
fMinIdx  = Result.glcDirect.WarmStartInfo.fMinIdx;

%if isempty(minPnts)
if ~feasible
   dBIG = 1E300;
   F = Result.glcDirect.WarmStartInfo.F(:);
   
   % F may not be a 0x1 matrix, but has to be a 0x0 matrix if empty
   if isempty(F)
     F = [];
   end
   C = Result.glcDirect.WarmStartInfo.C;
   G = Result.glcDirect.WarmStartInfo.G;

   if Inform < 90 
      Inform = 91; 
      ExitFlag = 7;
      Result.ExitText = [Result.ExitText ' No feasible point.'];
   end
   % Index sets for linear and nonlinear constraints
   % A   = Prob.A;        % Linear constraint matrix
   %b_L = Prob.b_L(:);   % Lower bounds, linear constraints
   %b_U = Prob.b_U(:);   % Upper bounds, linear constraints
   c_L = Prob.c_L(:);   % Lower bounds, nonlinear constraints
   c_U = Prob.c_U(:);   % Upper bounds, nonlinear constraints
   x_D = x_U-x_L;

   if isempty(b_L)
      leq     = [];
      b_L_idx = [];
      if isempty(b_U)
         b_U_idx = [];
      else
         b_U_idx = find(isfinite(b_U));
      end
   elseif isempty(b_U)
      leq     = [];
      b_U_idx = [];
      b_L_idx = find(isfinite(b_L));
   else
      leq     = b_L == b_U;
      b_L_idx = find(isfinite(b_L) & ~leq);
      b_U_idx = find(isfinite(b_U) & ~leq);
      leq     = find(leq);
   end
   dLin = size(Prob.A,1);
   if isempty(c_L)
      nleq    = [];
      c_L_idx = [];
      if isempty(c_U)
         c_U_idx = [];
      else
         c_U_idx = find(isfinite(c_U));
      end
   elseif isempty(c_U)
      nleq    = [];
      c_U_idx = [];
      c_L_idx = find(isfinite(c_L));
   else
      nleq    = c_L == c_U;
      c_L_idx = find(isfinite(c_L) & ~nleq);
      c_U_idx = find(isfinite(c_U) & ~nleq);
      nleq    = find(nleq);
   end

   LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
   NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
   mI    = LinI + NLinI;                     % Number of constraints
   LinE  = length(leq);
   NLinE = length(nleq);
   mE    = LinE + NLinE;
   m     = mI + mE;
   %g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
   %        b_L(leq);c_L(nleq)];
   g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
          cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
          bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
   N = size(G,2);
   if dLin > 0
      % Find point feasible w.r.t. linear constraints
      if LinI > 0
         if LinI > 1
            %ixI = all(G(1:LinI,:) - g_U(1:LinI)*ones(1,N) <= bTol);
            ixI = all(G(1:LinI,:) <= bTol);
         else
            %ixI = G(1:LinI,:) - g_U(1:LinI)*ones(1,N) <= bTol;
            ixI = G(1:LinI,:) <= bTol;
         end
         if LinE > 1
            %ixE = all(abs(G(mI+1:mI+LinE,:) - ...
            %          g_U(mI+1:mI+LinE)*ones(1,N)) <= bTol);
            ixE = all(abs(G(mI+1:mI+LinE,:)) <= bTol);
            ix  = ixE & ixI;
         elseif LinE == 1
            ixE = abs(G(mI+1:mI+LinE,:)) <= bTol;
            ix  = ixE & ixI;
         else
            ix  = ixI;
            ixE = [];
         end
      elseif LinE > 1
         ixE = all(abs(G(mI+1:mI+LinE,:)) <= bTol);
         ix  = ixE;
         ixI = [];
      elseif LinE == 1
         ixE = abs(G(mI+1:mI+LinE,:)) <= bTol;
         ix  = ixE;
         ixI = [];
      end
      if ~any(ix)
         if any(ixI) 
            ix = ixI;
            Afeas   = 2;
         elseif any(ixE) 
            ix = ixE;
            Afeas   = 2;
         else
            Afeas   = 0;
         end
      else
         Afeas   = 1;
      end
   else
      Afeas   = 3;
   end
   if mI + mE == 1
      if mI == 1
         Fpen = [G(1:mI,:) - g_T(1:mI)*ones(1,N)];
      else
         Fpen = [abs(G(mI+1:m,:))-g_T(mI+1:m)*ones(1,N)];
      end
   else
      if mI == 0
         Fpen = max(abs(G(mI+1:m,:))-g_T(mI+1:m)*ones(1,N));
      elseif mE == 0
         Fpen = max(G(1:mI,:) - g_T(1:mI)*ones(1,N));
      else
         Fpen = max([G(1:mI,:) - g_T(1:mI)*ones(1,N); ...
             abs(G(mI+1:m,:))-g_T(mI+1:m)*ones(1,N)]);
      end
   end
   if Afeas == 1
      ixx = find(ix);
      [fBest,iMin]=min(F(ixx)'+Fpen(ixx));
      fMinIdx = ixx(iMin);
   elseif Afeas == 2
      ixx = find(ix);
      [fBest,iMin]=min(Fpen(ixx));
      fMinIdx = ixx(iMin);
   elseif Afeas == 3
      [fBest,fMinIdx]=min(F'+Fpen);
   else
      [fBest,fMinIdx]=min(Fpen);
   end
   %Afeas
   %fMinIdx
   Result.minPnts = fMinIdx;
   xMin = tomsol(9,x_L,C(:,fMinIdx),x_D);    
   Result.x_k     = xMin;
   Result.f_k = F(fMinIdx);
   if Result.f_k ==  dBIG
      Result.f_k = nlp_f(xMin, Prob);  % Function value at xBest
   end

   if LinI + LinE > 0 & size(Result.x_k,2) == 1
      Result.Ax=Prob.A*Result.x_k; % Linear constraint value at all x_k;
   end
   Result.glcDirect.WarmStartInfo.fMinIdx = fMinIdx;
   %keyboard
elseif isempty(minPnts)
   % if minPnts [], no progress in this run, get best point from previous run
   F = Result.glcDirect.WarmStartInfo.F;
   % F may not be a 0x1 matrix, but has to be a 0x0 matrix if empty
   C = Result.glcDirect.WarmStartInfo.C;
   G = Result.glcDirect.WarmStartInfo.G;

   x_D            = x_U-x_L;
   Result.minPnts = fMinIdx;
   xMin           = tomsol(9,x_L,C(:,fMinIdx),x_D);    
   Result.x_k     = xMin;
   Result.f_k     = F(fMinIdx);
   Inform         = 10;
else
   Result.minPnts = minPnts;
   Result.x_k     = xMin;         % Best point(s) found
   Result.f_k     = glcfMin;      % Best function value
end
if length(Prob.c_L) + length(Prob.c_U)  > 0
   c_k=[];
   for i=1:size(xMin,2)
       cxx=nlp_c(xMin(:,i), Prob);
       c_k = [c_k  cxx(:)];
   end
   Result.c_k      = c_k; % Constraint value at x_k;
   Result.ConstrEv = Result.FuncEv;
end

Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

SaveWarmStartFile(Result.glcDirect.WarmStartInfo);

% Must set global n_f to see number of function evaluations in the GUI

global n_f
n_f = Result.FuncEv;

Result = endSolve(Prob,Result);


function SaveWarmStartFile(WarmStartInfo)

%     .nFunc
%     .nIter
%     .feasible
%     .fMinIdx
%     .s0
%     .glcfMin
%     .fMinEQ
%     .C
%     .F
%     .D
%     .G
%     .roc
%     .Tr
%     .ign
%     .iL
%     .iU
%     .Split
%     .t
%     .nFuncTot
%     .nRectTot
%     plus .Name

filename = 'glcDirectSave.mat';

nFunc    = WarmStartInfo.nFunc;
nIter    = WarmStartInfo.nIter;
feasible = WarmStartInfo.feasible;
fMinIdx  = WarmStartInfo.fMinIdx;
s0       = WarmStartInfo.s0;
glcfMin  = WarmStartInfo.glcfMin;
fMinEQ   = WarmStartInfo.fMinEQ;
C        = WarmStartInfo.C;
F        = WarmStartInfo.F;
D        = WarmStartInfo.D;
G        = WarmStartInfo.G;
roc      = WarmStartInfo.roc;
Tr       = WarmStartInfo.Tr;
ign      = WarmStartInfo.ign;
iL       = WarmStartInfo.iL;
iU       = WarmStartInfo.iU;
Split    = WarmStartInfo.Split;
t        = WarmStartInfo.t;
nFuncTot = WarmStartInfo.nFuncTot;
nRectTot = WarmStartInfo.nRectTot;
Name     = WarmStartInfo.Name;

try
   save(filename,'nFunc','nIter','feasible', ...
      'fMinIdx','s0','glcfMin','fMinEQ', ...
      'C', 'F','D','G','roc','Tr','ign', ...
      'iL','iU','Split','t','nFuncTot',...
      'nRectTot','Name');
catch
   warning('Failed to save warm start information to glcDirectSave.mat')
   disp(lasterr)
   disp('Warm start information is available in Result.glcDirect.WarmStartInfo')
end

% MODIFICATION LOG
%
% 050329 frhe File written, based on glcFast.m.
% 050330 frhe Some convFlag cases added. IntVars checks added.
% 050331 frhe Help modified.
% 050405 hkh  Move LinWeight, fEqual from Prob.GO to Prob.glcDirect
% 050405 hkh  Add fcALL to Prob.glcDirect, change AxFeas
% 050420 frhe ExitText messages changed
% 050420 frhe Call to mex changed and new parameters added
% 060327 hkh  Allow both 2 and 3 user function input arguments (3rd=varargin)
% 060327 hkh  Add output of Result.ConstrEv=Result.FuncEv, if mNonLin > 0
% 060814 med  FUNCS used for callbacks instead
% 061003 ango try-catch safe of save statement
% 061005 med  Minor help updates
% 070222 hkh  Revise IntVars handling, use new format, safe input to mex
% 070803 ango Mex call changed to use probPtr, not Prob
% 091002 hkh  Use WarmDefGLOBAL to generate warm start info in Prob structure
% 091002 hkh  Revise and correct comments about warm start
% 091002 hkh  Improve ExitText if WarmStart=1, similar to glbDirect
