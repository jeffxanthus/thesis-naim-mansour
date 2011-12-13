% glcSolve.m
%
% Solves general constrained mixed integer global optimization problems.
%
% glcSolve.m implements the algorithm DIRECT by Donald R. Jones presented
% in the paper "DIRECT", Encyclopedia of Optimization, Kluwer Academic
% Publishers, 2001.
% The algorithm is expanded to handle nonlinear and linear equalities, and
% linear inequalities.
%
% glcSolve solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Recommendation: Put the integers as the first variables.
% Put low range integers before large range integers.
% Linear constraints are specially treated.
% Equality constraints are added as penalties to the objective.
% Weights are computed automatically, assuming f(x) scaled to be roughly 1
% at optimum. Otherwise scale f(x).
%
% glcSolve will set f=100000 and not compute f(x) for any x where the linear
% constraints are not feasible. Any nonlinear constraints are computed.
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"
%
% Calling syntax:
%
% function Result = glcSolve(Prob, varargin )
%
% INPUT PARAMETERS
%
% Prob    Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say GLCF
%   FUNCS.c   The routine to compute the nonlinear constraint, say GLCC
%             A call to tomFiles could be used to set the names into
%             the Prob struct:
%             Prob = tomFiles(Prob,'GLCF',[],[],'GLCC');
%   x_L       Lower bounds for each element in x. Try to make a tight bound
%             Any lower bounds that are inf are changed to -10000
%   x_U       Upper bounds for each element in x. Try to make a tight bound
%             Any upper bounds that are inf are changed to  10000
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   PriLevOpt Print Level
%             0 = silent. 1 = some printing. 2 = print each iteration
%   WarmStart If true, >0, glcSolve reads the output from the last run
%             from the mat-file glcSave.mat, and continues from the last run.
%             NOTE: All rectangles that are fathomed in the previous run are
%             deleted. This saves space and computational time and enables
%             solving larger problems and more function evaluations to be done
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%
% ---------------------------------------
% glcDirect   Structure with DIRECT algorithm specific parameters. Fields used:
% ---------------------------------------
%
%   fcALL     =0 (Default). If linear constraints cannot be feasible anywhere
%             inside rectangle, skip f(x) and c(x) computation for middle point
%             =1 Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Do not update rates of change 
%             for the constraints.
%             =2 Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Update rates of change
%             constraints.
%
%   useRoC    =1 (Default). Use original Rate of Change (RoC) for constraints
%             to weight the constraint violations in selecting which rectangle
%             divide.
%             =0 Avoid RoC, giving equal weights to all constraint violations.
%             Suggested if difficulty to find feasible points. For problems
%             where linear constraints have been added among the nonlinear
%             (NOT RECOMMENDED; AVOID!!!), then option useRoc=0 has been 
%             successful, whereas useRoC completely fails. 
%             =2 Avoid RoC for linear constraints, giving weight one to these
%             constraint violations, whereas the nonlinear constraints use RoC.
%             =3 Use RoC for nonlinear constraints, but linear constraints are
%             not used to determine which rectangle to use.
%
%   BRANCH    =0 Divide rectangle by selecting the longest side, if ties use
%             the lowest index. This is the Jones DIRECT paper strategy.
%             =1 First branch the integer variables, selecting the variable
%             with the least splits. If all integer variables are split, split
%             on the continuous variables as in BRANCH=0.  DEFAULT!
%             Normally much more efficient than =0 for mixed-integer problems
%             =2 First branch the integer variables with 1,2 or 3 possible values,
%             e.g [0,1],[0,2] variables, selecting the variable with least splits. 
%             Then branch the other integer variables, selecting the variable
%             with the least splits. If all integer variables are split, split
%             on the continuous variables as in BRANCH=0. 
%             =3 Like =2, but use priorities on the variables, similar to
%             mipSolve, see Prob.MIP.VarWeight.
%   RECTIE    When minimizing the measure to find which new rectangle to try to 
%             get feasible, there are often ties, several rectangles have the same
%             minimum. RECTIE = 0 or 1 seems reasonable choices. Rectangles with low
%             index are often larger then the rectangles with higher index.
%             Selecting one of each type could help, but often =0 is fastest.
%             =0 Use the rectangle with value a, with lowest index (original)
%             =1 (Default): Use 1 of the smallest and 1 of largest rectangles
%             =2 Use the last rectangle with the same value a, not the 1st
%             =3 Use one of the smallest rectangles with same value a
%             =4 Use all rectangles with the same value a, not just the 1st
%
%   EqConFac  Weight factor for equality constraints when adding to objective
%             function f(x) (Default value 10). The weight is computed as
%             EqConFac/"right or left hand side constant value", e.g. if
%             the constraint is Ax <= b, the weight is EqConFac/b
%             If DIRECT just is pushing down the f(x) value instead of
%             fulfilling the equality constraints, increase EqConFac
%
%   AxFeas    Set nonzero to make glcSolve skip f(x) evaluations, when the
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
% ---------------------------------------
% optParam    Structure in Prob, Prob.optParam.
% ---------------------------------------
%             Defines optimization parameters. Fields used:
%   MaxIter   Maximal number of iterations, default max(5000,n*1000);
%   MaxFunc   Maximal number of function evaluations, default max(10000,n*2000)
%   IterPrint Print one line each iteration
%   cTol      Nonlinear constraint feasibility tolerance
%   bTol      Linear constraint feasibility tolerance
%   fGoal     Goal for function value, if empty not used
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%   eps_x     Convergence tolerance in x. All possible rectangles are
%             less than this tolerance (scaled to (0,1) )
%             See the output field maxTri.
%   EpsGlob   Global/local weight parameter, default 1E-4.
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
%           7 = Reached maxFunc or maxIter, NOT feasible
%           8 = Empty domain for integer variables
%  Inform   1 = Function value f is less than fGoal
%           2 = Absolute function value f is less than fTol, only if fGoal = 0
%            or Relative error in function value f is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           3 = Maximum number of iterations done
%           4 = Maximum number of function evaluations done
%           9 = Max CPU Time reached
%           91= Infeasible
%           99= Input error, see ExitFlag
%
% To make a warm start possible, glcSolve saves the following information in
% the file glcSave.mat (for internal solver use only):
%   C         Matrix with all rectangle centerpoints, in [0,1]-space.
%   D         Vector with distances from centerpoint to the vertices.
%   F         Vector with function values.
%   G         Matrix with constraint values for each point.
%   Name      Name of the problem. Used for security if doing warm start
%   Split     Split(i,j) = # splits along dimension i of rectangle j
%   T         T(i) is the number of times rectangle i has been trisected.
%   fMinEQ    sum(abs(infeasibilities)) for minimum points, 0 if no equalities
%   fMinIdx   Indices of the currently best points
%   feasible  Flag indicating if a feasible point has been found.
%   glcfMin   Best function value found at a feasible point.
%   iL        iL(i,j) is the lower bound for rect. j in integer dim. I(i)
%   iU        iU(i,j) is the upper bound for rect. j in integer dim. I(i)
%   ignoreIdx Rectangles to be ignored in the rect. selection proceedure.
%   s         s(j) is the sum of observed rates of change for constraint j.
%   s_0       s_0 is used as s(0)
%   t         t(i) is the total # splits along dimension i.
%   SubRes    Additional output from nlp_f, if suboptimization done
%
% USAGE:
%
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
%      The name of the problem is set to "GLCF Test"
%
%      Prob   = glcAssign('GLCF',x_L,x_U,'GLCF Test',A,b_L,b_U,'GLCC',c_L,c_U);
%      Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%      % Default values are now set for PriLevOpt, and structure optParam
%      % To change a value, an example is shown on the next line
%      Prob.optParam.MaxFunc = 500; % Change max number of function evaluations
%
%      If there are integer variables, they may be set as additional input
%      to glcAssign, or directly as the next line shows:
%      Prob.MIP.IntVars = [1 3];  % 1st and third variables are integers
%
%      Result = tomRun('glcSolve',Prob,2); % Print level 2
%
% The user function GLCF is written as
%
%      function f = GLCF(x, Prob)
%      u = Prob.user.u; W = Prob.user.W;
%      f = "some function of x, u and W"
%
% It is also possible to use the function format
%      function f = GLCF(x)
% or sending the extra parameters as additional input
%
%      Result = glcSolve(Prob,u,W);
%
%      function f = GLCF(x,Prob,u,W)
%
% NOTE! If additional parameters are sent, Prob must be the second input
% parameter to GLCF (and GLCC)
%
% The user function GLCC, computing the nonlinear constraints, is written as
%
%      function c = GLCC(x, Prob)
%      u = Prob.user.u; W = Prob.user.W;
%      c = "some vector function of x, V and W"
%
% To make a restart, just set the restart flag, and call glcSolve once again:
%
%      Prob.WarmStart = 1;
%      Result = tomRun('glcSolve',Prob,2); % Print level 2
%
%      % It is recommended to do frequent warm starts if to do many
%      % function evaluations
%      Prob.optParam.MaxFunc = 10000;
%      Result = tomRun('glcSolve',Prob,2); % 1st call without warm start
%      Set a high the maximum number of iterations
%      Prob.optParam.MaxIter = 100000; % will never be reached
%      Prob.WarmStart        = 1;
%      for i=1:9
%          Result = tomRun('glcSolve',Prob,2); % Print level 2
%      end
%      Now around 100000 rectangles have been tried
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
% upper bounds, as much as possible.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Feb 15, 1999.   Last modified Feb 22, 2007.

function Result = glcSolve(Prob, varargin )

if nargin < 1
   error('glcSolve needs input structure Prob');
end

solvType=checkType('glc');

Prob=ProbCheck(Prob,'glcSolve',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = Prob.MaxCPU;

% Pick up input parameters from the Prob structure:
x_L = Prob.x_L(:);   % Lower bounds
x_U = Prob.x_U(:);   % Upper bounds

if isempty(x_L) | isempty(x_U) | any(isinf(x_L)) | any(isinf(x_U))
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
      error('glcSolve: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

PriLev     = Prob.PriLevOpt;          % Print level
LargeScale = Prob.LargeScale;

MaxFunc    = Prob.optParam.MaxFunc;   % Maximum number of function evaluations
MaxIter    = Prob.optParam.MaxIter;   % Number of iterations
EpsGlob    = Prob.optParam.EpsGlob;   % global/local weight parameter. 
cTol       = Prob.optParam.cTol;      % Constraint feasibility tolerance
bTol       = Prob.optParam.bTol;      % Linear Constraint feasibility tolerance
IterPrint  = Prob.optParam.IterPrint; % Print short information each iteration
fGoal      = Prob.optParam.fGoal;     % Goal for f(x).
fTol       = Prob.optParam.eps_f;     % Relative tolerance for fGoal
xTol       = Prob.optParam.eps_x;     % Tolerance for rectangle sizes

% Safeguard
if isempty(MaxIter) | MaxIter < 0
   MaxIter = max(5000,n*1000);
end
if isempty(MaxFunc) | MaxFunc < 0
   MaxFunc = max(10000,n*2000);
end
                                     % (scaled to (0,1) )
glcfOld = Inf;
if isempty(LargeScale) 
   if MaxFunc >= 10000
      LargeScale = 1; 
   else
      LargeScale = 0; 
   end
end

% glcDirect parameters
GD = DefPar(Prob,'glcDirect',[]);

fcALL     = DefPar(GD,'fcALL',0); 
useRoC    = DefPar(GD,'useRoC',1); 
BRANCH    = DefPar(GD,'BRANCH',1); 
RECTIE    = DefPar(GD,'RECTIE',1); 
EqConFac  = DefPar(GD,'EqConFac',10); 
FSTAR     = DefPar(GD,'FSTAR ',0); 
AxFeas    = DefPar(GD,'AxFeas',0); 
fEqual    = DefPar(GD,'fEqual',1E-10); 
LinWeight = DefPar(GD,'LinWeight',0.1); 
alfa      = max(0,min(1,DefPar(GD,'alpha',0.9))); 
AvIter    = max(1,DefPar(GD,'AvIter',50)); 

AxFeas    = min(fcALL,AxFeas);  % AxFeas must be 0 if fcALL is 0

%useRoC = 1;
%EqConFac  = 10;
%EqConFac  = 0.0001;
%EqConFac  = 100000;
%BRANCH = 3;
%RECTIE = 1;
%FSTAR  = 0;
if PriLev > 0 | IterPrint > 0
   fprintf('useRoC %d EqConFac %d BRANCH %d RECTIE %d FSTAR %d\n',...
            useRoC,EqConFac,BRANCH,RECTIE,FSTAR);
end

SaveRes = DefPar(Prob,'SaveRes',0);
SubRes  = [];
Res     = [];


nFunc    = 0;        % Function evaluation counter
nFuncTot = 0;        % Total Function evaluation counter
nRectTot = 0;        % Total number of rectangles tried

Result                 = ResultDef(Prob);
Result.Solver          = 'glcSolve';
Result.SolverAlgorithm = 'DIRECT - Lipschitzian Optimization';


A   = Prob.A;        % Linear constraint matrix
b_L = Prob.b_L(:);   % Lower bounds, linear constraints
b_U = Prob.b_U(:);   % Upper bounds, linear constraints
c_L = Prob.c_L(:);   % Lower bounds, nonlinear constraints
c_U = Prob.c_U(:);   % Upper bounds, nonlinear constraints


isMIP   = ~isempty(IntVars); % isMIP is true if there are integer variables
I_logic = zeros(n,1);        % Logic set for integer/continuous variables
I_logic(IntVars) = 1;
if isMIP & BRANCH == 3
   VarWeight = DefPar(Prob.MIP,'VarWeight',[]); 
   if isempty(VarWeight)
      VarWeight = ones(length(IntVars),1);
   end
end

% Index sets for linear and nonlinear constraints
beta = [];
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
   if ~isempty(leq)
      if EqConFac == 1
         beta = ones(length(leq),1);
      else
         v       = b_L(leq);
         v(v==0) = 1;
         % Weight on infeasibility addition
         beta    = abs(EqConFac./v);
      end
   end
end
dLin = size(Prob.A,1);
%dLin = 0;
dSum = 0;
dBIG = 100000;
if dLin > 0
   d1 = length(Prob.b_L);
   % Must adjust bounds to have equal length
   if d1 < dLin
      b_L = [b_L;-inf*ones(dLin-d1,1)];
   end
   d1 = length(Prob.b_U);
   if d1 < dLin
      b_U = [b_U;inf*ones(dLin-d1,1)];
   end
end
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
   if ~isempty(nleq)
      if EqConFac == 1
         beta = [beta;ones(length(nleq),1)];
      else
         v       = c_L(nleq);
         v(v==0) = 1;
         % Weight on infeasibility addition
         beta    = [beta;abs(EqConFac./v)];
      end
   end
end
fEqualW = sum(beta)*fEqual; % Balance equalities due to beta weight


%b_L_idx = find(~isempty(b_L) & isfinite(b_L));
%b_U_idx = find(~isempty(b_U) & isfinite(b_U));
%c_L_idx = find(~isempty(c_L) & isfinite(c_L));
%c_U_idx = find(~isempty(c_U) & isfinite(c_U));

% Check if there is a chance to find a feasible point
if any(x_L > x_U)
   fprintf('\n\n');
   fprintf('Error in glcSolve, upper variable bounds below lower bounds\n');
   Result.ExitFlag = 2;
   Result.Inform   = 99;
   Result.ExitText =str2mat('glcSolve solves box-bounded problems' ...
            ,'Found some upper bound below lower bounds');
   Result.maxTri = [];
   Result=endSolve(Prob,Result);
   return;
end
x_D   = x_U-x_L;
nFix  = sum(x_D == 0);
nCont = sum(x_D > 0 & ~I_logic);

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
      fIP = nlp_f(xIP, Prob, varargin{:});  % Function value at x
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
   if ~isempty(A)
      Ax = A*xIP;
      if any(xIP < x_L | xIP > x_U) | ...
         any(b_L-Ax > bTol) | any(Ax-b_U > bTol)
         if PriLev > 0 | IterPrint > 0
            disp('Input Prob.MIP.xIP is not linear feasible');
            disp('Reject value of Prob.MIP.fIP');
         end
         fIP = Inf;     
         % xIP = [];     
      end
   end
end

%
%  STEP 1, INITIALIZATION
%

if Prob.WarmStart
   % Restart with values from previous run.

   load('glcSave.mat','Name','C','F','T','D','G','iL','iU','s_0','s','t',...
         'ignoreIdx','feasible','Split','glcfMin','fMinIdx','fMinEQ',...
         'SubRes','nFuncTot','nRectTot');
   Name1 = Prob.Name;               % Name for the problem
   if strcmp(Name1,Name)
      LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
      NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
      mI    = LinI + NLinI;                     % Number of constraints
      LinE  = length(leq);
      NLinE = length(nleq);
      mE    = LinE + NLinE;
      m     = mI + mE;
      nFuncOld = size(C,2);
      nOld     = nFuncOld;
      maxFDim  = MaxFunc + max(100,10*n);   
      if isempty(ignoreIdx)
         F        = [F,zeros(1,maxFDim)]; 
         D        = [D,zeros(1,maxFDim)]; 
         if LargeScale
            T     = [T;sparse(maxFDim,1)]; 
            C     = [C,sparse(n,maxFDim)];
            Split = [Split,sparse(n,maxFDim)]; 
            if isMIP
               IntVarsLen = length(IntVars);
               iL         = [iL,sparse(IntVarsLen,maxFDim)];
               iU         = [iU,sparse(IntVarsLen,maxFDim)];
            end
         else
            % T(i), number of times rectangle i has been trisected
            %T = [T;zeros(nFuncOld+MaxFunc-length(T),1)]; 
            T     = [T;zeros(maxFDim,1)]; 
            C     = [C,zeros(n,maxFDim)];
            Split = [Split,zeros(n,maxFDim)]; 
            if isMIP
               IntVarsLen = length(IntVars);
               iL         = [iL,zeros(IntVarsLen,maxFDim)];
               iU         = [iU,zeros(IntVarsLen,maxFDim)];
            end
         end
      else
         ix            = ones(nFuncOld,1);
         ix(ignoreIdx) = 0;
         if isnan(glcfMin)
            ix            = find(ix);
            fMinIdx       = 1;
         else
            ix(fMinIdx)   = 1;
            ix            = find(ix);
            for i = 1:length(fMinIdx)
                fMinIdx(i) = find(fMinIdx(i)==ix);
            end
         end
         F        = [F(ix),zeros(1,maxFDim)]; 
         D        = [D(ix),zeros(1,maxFDim)]; 
         if LargeScale
            T     = [T(ix);sparse(maxFDim,1)]; 
            C     = [C(:,ix),sparse(n,maxFDim)];
            Split = [Split(:,ix),sparse(n,maxFDim)]; 
            if isMIP
               IntVarsLen = length(IntVars);
               iL         = [iL(:,ix),sparse(IntVarsLen,maxFDim)];
               iU         = [iU(:,ix),sparse(IntVarsLen,maxFDim)];
            end
         else
            % T(i), number of times rectangle i has been trisected
            %T = [T(ix);zeros(nFuncOld+MaxFunc-length(T),1)]; 
            T     = [T(ix);zeros(maxFDim,1)]; 
            C     = [C(:,ix),zeros(n,maxFDim)];
            Split = [Split(:,ix),zeros(n,maxFDim)]; 
            if isMIP
               IntVarsLen = length(IntVars);
               iL         = [iL(:,ix),zeros(IntVarsLen,maxFDim)];
               iU         = [iU(:,ix),zeros(IntVarsLen,maxFDim)];
            end
         end
         nFuncOld  = length(ix); % Now fathomized rectangles are deleted
      end
      if m > 0
         % Place all constraints in one vector gx with upper bounds g_U
         g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
               b_L(leq);c_L(nleq)];
         g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
                cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
                bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
         if isempty(ignoreIdx)
            if LargeScale
               G          = [G,sparse(m,maxFDim)];
            else
               G          = [G,zeros(m,maxFDim)];
            end
         else
            if LargeScale
               G          = [G(:,ix),sparse(m,maxFDim)];
            else
               G          = [G(:,ix),zeros(m,maxFDim)];
            end
         end
      else
         rE = 0;
      end
      ignoreIdx = [];
      AvFunc  = nFuncTot;
      nFunc   = nFuncOld;            % RESET nFunc
      MaxFunc = MaxFunc + nFuncOld;  % RESET MaxFunc
  
      if PriLev > 0 | IterPrint > 0
         fprintf('\n');
         fprintf('Restart with %d sampled points from previous run',nFuncOld);
         fprintf(', %d ignored',nOld-nFuncOld);
         fprintf('\n');
      end
      if ~isinf(glcfMin) 
         % All points i with F(i)=glcfMin, and feasible
         xBest = tomsol(9,x_L,full(C(:,fMinIdx)),x_D);    
      else
         xBest = [];
      end
   else
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files glcSave.mat?\n');
      end
   end
end % if WarmStart

if ~Prob.WarmStart
   Name = deblank(Prob.Name);  % Problem name
   % No restart, set first point to center of the unit hypercube.
   nFuncOld = 0;
   maxFDim = MaxFunc + max(100,10*n);   

   % SAMPLE THE CENTERPOINT OF THE ENTIRE SPACE.
   if LargeScale
      C = sparse(n,maxFDim);
   else
      C = zeros(n,maxFDim);
   end
   C(:,1) = ones(n,1)./2;     % Matrix with all rectangle centerpoints.
   F      = zeros(1,MaxFunc); % Vector with function values
   % All C_coordinates refers to the n-dimensional hypercube. 
   % Split(i,j) is the number of times rectangle j has  been split along dimension i. Split(i,j) is set to
   % Inf if rect j has sidelength 0 for integer variable i. 
   if LargeScale
      Split = sparse(n,maxFDim); 
      T     = sparse(MaxFunc,1); % T(i), number of times rectangle i has been trisected
   else
      %Split = zeros(n,1); 
      Split = zeros(n,maxFDim); 
      T     = zeros(MaxFunc,1); % T(i), number of times rectangle i has been trisected
   end

   % IF ALL VARIABLES ARE INTEGERS, THEN IT IS POSSIBLE FOR A RECTANGLE
   % TO BE REDUCED TO A SINGLE POINT. IF THIS HAPPENS, THE RECTANGLE SHOULD
   % BE IGNORED IN THE RECTANGLE SELECTION PROCEEDURE.
   ignoreIdx = []; 

   % Fixed variables should not be considered for the splitting procedure
   idx = find(x_L==x_U);
   % Eliminates risk choosing these dims as splitting dims for any rectangle 
   if length(idx) > 0
      Split(idx,1) = Inf; 
      if length(idx)==n
         ignoreIdx = [ignoreIdx;1];
      end
   end
   % HKH Move nFunc definition to first
   nFunc = nFunc+1;          % Number of function values

   if isMIP
      % iL(i,j) is the lower bound for rectangle j in integer dimension I(i).
      % iU(i,j) is the upper bound for rectangle j in integer dimension I(i).
      IntVarsLen = length(IntVars);
      if LargeScale
         iL = sparse(IntVarsLen,maxFDim);
         iU = sparse(IntVarsLen,maxFDim);
      else
         iL = zeros(IntVarsLen,maxFDim);
         iU = zeros(IntVarsLen,maxFDim);
      end
      iL(:,1) =  ceil(x_L(IntVars)); 
      iU(:,1) = floor(x_U(IntVars)); 

      IUL = full(iU(:,1)-iL(:,1));
      if any(IUL < 0)
         disp('Error in glcSolve, empty domain for integer variables:')
         tmpidx = find(IUL<0);
         IUL(tmpidx)
         Result.ExitFlag = 8;
         Result.ExitText = 'Empty domain for integer variables';
         Result.maxTri   = [];
         Result.Inform   = 99;
         Result=endSolve(Prob,Result);
         return;
      end
      x_mid  = full(floor( (iU(:,1)+iL(:,1))/2 ));
      %x_mid  = floor( (x_U(I)+x_L(I))/2 );
      C(IntVars,1) = (x_mid-x_L(IntVars))./(max(x_D(IntVars),1E-20));
      tmpidx    = find(IUL==0);
      % Eliminate risk choosing these dims as splitting dim for any rectangle
      if ~isempty(tmpidx)
         Split(IntVars(tmpidx),1) = Inf; 
         if length(tmpidx)==n
            ignoreIdx = [ignoreIdx;1];
         end 
      end
   else
      iL = [];
      iU = [];
   end

   % Transform C to original search space
   %x     = x_L + C.*x_D;  
   x    = tomsol(9,x_L, full(C(:,1)),x_D); 
   if dLin > 0
      Ax = A*x;
      %Afeas = sum(max(0,max(b_L-bTol-Ax,Ax-bTol-b_U)));
      %format long
      %[b_L, Ax, b_U]
   else
      Ax=[];
      %Afeas = 0;
   end   
   % if Afeas == 0
   if 1  % Always compute the first point in space
      if SaveRes
         Prob.P = nFunc;
         [f,Res]  = nlp_f(x, Prob, varargin{:});  % Save additional result
      else
         f        = nlp_f(x, Prob, varargin{:});  % Function value at x
      end
   else
      f     = dBIG;
      dSum  = dSum +1;
   end
   D     = zeros(1,maxFDim); 
   D(1)  = sqrt(n-nFix)/2;   % Vector with distances from centerpoint to the vertices

   % IF THE CENTER IS FEASIBLE, SET x_min EQUAL TO THE CENTERPOINT AND
   % glcfMin EQUAL TO THE OBJECTIVE FUNCTION VALUE AT THIS POINT.

   LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
   NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
   mI    = LinI + NLinI;                     % Number of constraints
   LinE  = length(leq);
   NLinE = length(nleq);
   mE    = LinE + NLinE;
   m     = mI + mE;
   if NLinI+NLinE > 0 
      cx = nlp_c(x, Prob, varargin{:});
      cx = cx(:);
   else
      cx=[];
   end
   %if LinE > 0
   %   % Use penalty approach for linear equality constraints
   %   %f = f + sum(abs(Ax(leq,:)-b_L(leq)))/bTol;
   %   %f = f + beta*sum(abs(Ax(leq,:)-b_L(leq)));
   %end
   %if NLinE > 0
   %   % Use penalty approach for nonlinear equality constraints
   %   %f = f + sum(abs(cx(nleq,:)-c_L(nleq)))/cTol;
   %   %f = f + beta*sum(abs(cx(nleq,:)-c_L(nleq)));
   %end

   F(1,1) = f;  % Vector with function values

   if m > 0
      % Place all constraints in one vector gx with upper bounds g_U
      g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
              b_L(leq);c_L(nleq)];
      g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
             cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
             bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
      gx  = [-Ax(b_L_idx);Ax(b_U_idx);-cx(c_L_idx);cx(c_U_idx); ...
              Ax(leq);cx(nleq)]; 
      tmp = gx - g_U;
      if all( tmp(1:mI) <= g_T(1:mI) ) % Initial point is feasible ?
         feasible = 1; % feasible point has been found
      else
         feasible = 0;
      end
      if LargeScale
         G          = sparse(m,maxFDim);
      else
         G          = zeros(m,maxFDim);
      end
      %G      = zeros(m,MaxFunc + 20);
      G(:,1) = tmp; % Subtract g_U to get g(x)<=0 (or g(x) = 0 for ineq.)
      % G is a matrix with constraint values for each point.
      % G(i,j) is the value of constraint i at point j.
   else
      G = [];
      feasible = 1;
   end
   fMinIdx = 1;
   xBest   = x;
   SubRes  = Res;
   if mE > 0
      fMinEQ = full(sum(beta.*abs(G(mI+1:m,nFunc))));
   else
      fMinEQ = 0;
   end
   if feasible
      glcfMin = f+fMinEQ;
   else
      glcfMin = NaN;
   end

   % SET s(j)=0 FOR j=1,2,...,m.
   if m > 0
      s_0 = 0;          % Used as s(0).
      % s(j) is the sum of observed rates of change for constraint j.
      s   = zeros(m,1); 
      if LinI + LinE > 0
         if n == 1
            z = sqrt(A.^2');
         else
            z = sqrt(sum(A.^2'))';
         end
      end
      if LinI > 0
         s(1:LinI) = LinWeight*[z(b_L_idx);z(b_U_idx)];
      end
      if LinE > 0
         s(mI+1:mI+LinE) = LinWeight*z(leq);
      end
   else
      s_0 = [];
      s = [];
      rE = 0;
   end

   AvFunc = 0;

   % SET t(i)=0 FOR i=1,2,...,n.
   % t(i) is the number of times a rectangle has been split along dimension i.
   t = zeros(n,1); 
end % if not restart

Iter      = 0; % Iteration counter
convFlag  = 0;
cpumax    = 0;
TIME0     = Prob.TIME0;
NowFeas   = 0; % Avoid crash if IterPrint == 1 | PriLev >0, and m == 0 
if BRANCH == 2
   if isempty(IntVars)
      IX1 = [];
      IX2 = [];
   else
      IX1 = find((iU(:,1)-iL(:,1))<=2);
      IX2 = find((iU(:,1)-iL(:,1))> 2);
   end
end
   

while Iter < MaxIter & nFunc < MaxFunc & convFlag == 0 

   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   Iter = Iter+1;

   %
   %  STEP 2, SELECT RECTANGLES 
   %

   if length(ignoreIdx)==nFunc % If all rectangles are fathomed
      break
   end

   % COMPUTE c(j) VALUES USING THE CURRENT VALUES OF s_0 AND s(j), j=1,2,...,m  
   %if m > 0
   %   c = s_0./(max(s,1E-5));
   %   cGD = (c(1:mI)'*max(G(1:mI,1:nFunc)-cTol,0))./D;
   %   %cGD = max(c*ones(1,nFunc).*max(G(:,1:nFunc)-...
   %   %          cTol,0))./D;
   %end

   S = []; % The set of rectangles selected for trisection

   ixD = find(D(1:nFunc) == 0);
   D(ixD) = 1;
   if ~feasible
      % Weighted sum of violations in all inequality constraints
      if useRoC == 1
         c   = s_0./(max(s,1E-5));
         cGD = (c(1:mI)'*max(G(1:mI,1:nFunc)-cTol,0))./D(1:nFunc);
      elseif useRoC == 0
         % Option 1: Avoid rate of change weighting
         cGD = sum(max(G(1:mI,1:nFunc)-cTol,0),1)./D(1:nFunc);
      elseif useRoC == 2
         % No weight on linear constraints. RoC on nonlinear constraints
         c   = s_0./(max(s,1E-5));
         cGD = (sum(max(G(1:LinI,1:nFunc)-cTol,0),1) + ...
               c(LinI+1:mI)'*max(G(LinI+1:mI,1:nFunc)-cTol,0))./D(1:nFunc);
      else   % useRoC == 3
         % If any nonlinear, skip linear constraints. Use RoC on nonlinear constraints
         % If no nonlinear, use linear constraints without RoC
         c   = s_0./(max(s,1E-5));
         if mI > LinI
            cGD = c(LinI+1:mI)'*max(G(LinI+1:mI,1:nFunc)-cTol,0)./D(1:nFunc);
         else
            cGD = sum(max(G(1:LinI,1:nFunc)-cTol,0),1)./D(1:nFunc);
         end
      end
%if Iter <= 10 
%   for ii=1:nFunc
%       if ~any(ii==ignoreIdx)
%          fprintf('%d %f x ',ii,cGD(ii))
%          fprintf('%f ',C(:,ii))
%          fprintf(' G ')
%          fprintf('%f ',G(:,ii))
%          %fprintf(' D ')
%          %fprintf('%f ',D(ii))
%          fprintf(' c ')
%          fprintf('%f ',c)
%          fprintf('\n')
%       end
%   end
%end
      % IF A FEASIBLE TRIANGLE HAS NOT BEEN FOUND, SELECT THE RECTANGLE THAT
      % MINIMIZES THE RATE OF CHANGE REQUIRED TO BRING THE WEIGHTED CONSTRAINT
      % VIOLATIONS TO ZERO.
      cGD(ixD) = Inf;

      if ~isempty(ignoreIdx)
         ix0   =ones(length(cGD),1);
         ix0(ignoreIdx)=0;
         ix    = find(ix0);
         [a b] = min( cGD(ix) );
         if ~isfinite(a)
            disp(' No feasible integer point exist')
            S = [];
            break;
         end
         if RECTIE == 0
            % Basic: Use the rectangle with value a, with lowest index
            S = ix(b); 
         elseif RECTIE == 1
            % Option 1: Use one of the smallest and one of the largest rectangles
            b = find(cGD(ix)==a);
            if length(b) > 1
               [a ix0] = sort(D(ix(b)));
               S = ix([b(ix0(1)),b(ix0(end))]);
            else
               S = ix(b); 
            end
         elseif RECTIE == 2
            % Option 2: Use the last rectangle with same value a, not the 1st
            b = find(cGD(ix)==a);
            S = ix(b(end));
         elseif RECTIE == 3
            % Option 3: Use one of the smallest rectangles with same value a
            b = find(cGD(ix)==a);
            [a ix0] = min(D(ix(b)));
            S = ix(b(ix0(end)));
         elseif RECTIE == 4
            % Option 4: Use all rectangles with same value a, not just the 1st
            b = find(cGD(ix)==a);
            S = ix(b); 
         end
      else
         % Empty ignoreIdx, no rectangles deleted
         [a S] = min( cGD );
         if RECTIE == 0
            % Basic: Use the rectangle with value a, with lowest index
            % This is S
         elseif RECTIE == 1
            % Option 1: Use one of the smallest and one of the largest rectangles
            b = find(cGD==a);
            if length(b) > 1
               [a ix0] = sort(D(b));
               S = [b(ix0(1)),b(ix0(end))];
            else
               S = b; 
            end
         elseif RECTIE == 2
            % Option 2: Use the last rectangle with the same value a, not the 1st
            b = find(cGD==a);
            S = b(end);
         elseif RECTIE == 3
            % Option 3: Use one of the smallest rectangles with same value a
            b = find(cGD==a);
            [a ix0] = min(D(b));
            S = b(ix0(end));
         elseif RECTIE == 4
            % Option 4: Use all rectangles with the same value a, not just the 1st
            S = find(cGD==a);
         end
      end
%     if m > 10000
%     %if NLinI > 1
%        cGD = max(c*ones(1,nFunc).*max(G(:,1:nFunc)-...
%                  cTol,0))./D(1:nFunc);
%        if ~isempty(ignoreIdx)
%           [a b] = min( cGD(ix) );
%           b = ix(b);
%           if S~=b
%              S = [S,b];
%              S
%           end
%        else
%           [a b] = min( cGD );
%           if S~=b
%              S = [S,b];
%              S
%           end
%        end
%     end
%     if NLinI > 10000  % Do not use this idea
%        % Add extra rectangles when having several nonlinear constraints
%        for i = 1:NLinI
%            j = LinI + i;
%            cGD = max(G(j,1:nFunc),0)./D(1:nFunc);
%            iy = cGD > 0;
%            if ~isempty(ignoreIdx)
%               iy(ignoreIdx) = 0;
%            end
%            iy = find(iy);
%            if ~isempty(iy)
%               [a b] = min( cGD(iy) );
%               b = iy(b);
%               if ~any(b==S)
%                  S = [S,b];
%               else
%                  b = find( abs(cGD(iy) - a) < 1E-12 );
%                  if length(b) == 1, break; end
%                  b = iy(b);
%                  for k = 1:length(b)
%                      if ~any(b(k)==S)
%                         S = [S,b(k)];
%                         break;
%                      end
%                  end
%               end
%            end
%        end
%        S
%     end

      % [a b] = min( cGD );
      % % How to treat the integer variables????
      % if ~isempty(ignoreIdx)
      %    cGDtmp = cGD;
      %    % rectangle b is fathomed i.e. b is a single point and can't be split
      %    while ismember(b,ignoreIdx) 
      %       cGDtmp(b) = Inf;
      %       [a b] = min( cGDtmp );
      %       if ~isfinite(a)
      %          disp(' No feasible integer point exist')
      %          S = [];
      %          break;
      %       end
      %    end
      %    S = b;
      % else
      %    S = b;
      % end   
   else
      % ON THE OTHER HAND, IF A FEASIBLE POINT HAS BEEN FOUND, IDENTIFY THE
      % SET OF RECTANGLES THAT PARTICIPATE IN THE LOWER ENVELOPE. LET S BE
      % THE SET OF SELECTED RECTANGLES.
      if m > 0
         c = s_0./(max(s,1E-5));
         % Feasible w.r.t. inequalties - Include equalities into cGD
    
         if useRoC == 1
            % Use RoC on both linear and nonlinear constraints
            if mI == 0
               cGD = (c'*abs(G(:,1:nFunc)))./D(1:nFunc); 
            elseif mE == 0
               cGD = (c'*max(G(:,1:nFunc)-cTol,0))./D(1:nFunc); 
            else
               cGD = (c(1:mI)'*max(G(1:mI,1:nFunc)-cTol,0) + ...
                      c(mI+1:m)'*abs(G(mI+1:m,1:nFunc)))./D(1:nFunc); 
            end
         elseif useRoC == 0
            % No weights on linear and nonlinear constraints
            if mI == 0
               cGD = sum(abs(G(:,1:nFunc)),1)./D(1:nFunc); 
            elseif mE == 0
               cGD = sum(max(G(:,1:nFunc)-cTol,0),1)./D(1:nFunc); 
            else
               cGD = (sum(max(G(1:mI,1:nFunc)-cTol,0),1) + ...
                     sum(abs(G(mI+1:m,1:nFunc)),1))./D(1:nFunc); 
            end
         elseif useRoC == 2
            % No weight on linear constraints. RoC on nonlinear constraints
            if mI == 0
               cGD = (sum(abs(G(1:LinE,1:nFunc)),1) + ...
                      c(LinE+1:m)'*abs(G(LinE+1:m,1:nFunc)))./D(1:nFunc); 
            elseif mE == 0
               cGD = (sum(max(G(1:LinI,1:nFunc)-cTol,0),1) + ...
                      c(LinI+1:m)'*max(G(LinI+1:m,1:nFunc)-cTol,0))./D(1:nFunc); 
            else
               cGD = (sum(abs(G(mI+1:mI+LinE,1:nFunc)),1) + ...
                      c(mI+LinE+1:m)'*abs(G(mI+LinE+1:m,1:nFunc)) + ...
                      sum(max(G(1:LinI,1:nFunc)-cTol,0),1) + ...
                      c(LinI+1:mI)'*max(G(LinI+1:mI,1:nFunc)-cTol,0))./D(1:nFunc); 
            end
         else
            % If any nonlinear, skip linear constraints. Use RoC on nonlinear constraints
            % If no nonlinear, use linear constraints without RoC
            if mI == 0 & LinE == mE
               cGD = sum(abs(G(1:LinE,1:nFunc)),1)./D(1:nFunc); 
            elseif mI == 0 
               cGD = c(LinE+1:m)'*abs(G(LinE+1:m,1:nFunc))./D(1:nFunc); 
            elseif mE == 0 & LinI == mI
               cGD = sum(max(G(1:LinI,1:nFunc)-cTol,0),1)./D(1:nFunc); 
            elseif mE == 0
               cGD = c(LinI+1:m)'*max(G(LinI+1:m,1:nFunc)-cTol,0)./D(1:nFunc); 
            elseif LinI == mI
               cGD = c(mI+LinE+1:m)'*abs(G(mI+LinE+1:m,1:nFunc))./D(1:nFunc); 
            elseif LinE == mE
               cGD = c(LinI+1:mI)'*max(G(LinI+1:mI,1:nFunc)-cTol,0)./D(1:nFunc); 
            else
               cGD = (c(mI+LinE+1:m)'*abs(G(mI+LinE+1:m,1:nFunc)) + ...
                      c(LinI+1:mI)'*max(G(LinI+1:mI,1:nFunc)-cTol,0))./D(1:nFunc); 
            end
         end
%if Iter == 8
%   for ii=1:nFunc
%       if ~any(ii==ignoreIdx)
%          fprintf('%d %f x ',ii,cGD(ii))
%          fprintf('%f ',C(:,ii))
%          fprintf(' G ')
%          fprintf('%f ',G(:,ii))
%          %fprintf(' D ')
%          %fprintf('%f ',D(ii))
%          fprintf(' c ')
%          fprintf('%f ',c)
%          fprintf('\n')
%       end
%   end
%end
         %cGD = max(c*ones(1,nFunc).*max(G(:,1:nFunc)-...
         %          cTol,0))./D;
      end
      Epsilon = max(EpsGlob*abs(glcfMin),1E-8);

      if FSTAR == 0
         % f_star = glcfMin-Epsilon;
         % NEW !!! 
         f_star = min(fIP,glcfMin-Epsilon);
      elseif FSTAR == 1
         f_star = min(fIP,glcfMin-fMinEQ(1)-Epsilon);
      elseif FSTAR == 2
         if mE > 0 & fMinEQ > fEqualW
            f_star = min(fIP,min(F)-Epsilon);
         else
            f_star = min(fIP,glcfMin-fMinEQ(1)-Epsilon);
         end
      end

      %r_previous=Inf;
      r_previous=[];
      idx1=[];
      while 1
         f_star_max = -Inf;
         if m > 0
            h      = max((F(1:nFunc)-f_star),0)./D(1:nFunc) + cGD;
         else
            h      = max((F(1:nFunc)-f_star),0)./D(1:nFunc);
         end
         h(ixD) = Inf;

         tmpset1=ones(1,nFunc);
         tmpset1(idx1)      = 0; % not consider rects put in S in previous turn
         tmpset1(ignoreIdx) = 0; % not consider rects being fathomed
         tmpidx=find(tmpset1);
 
         if ~isempty(tmpidx)
            h_min  = min(h(tmpidx));
            idx1   = tmpidx(find(h(tmpidx)==h_min));
         else
            idx1       = [];
         end

         if length(idx1) > 1
            if any(f_star > F(idx1))
               % Choose rectangle being constant for the lowest value of f_star
               [a b] = max(f_star-F(idx1));
               r = idx1(b);
            else
               % Choose the rectangle with the flatest slope.
               [a b] = max(D(idx1));
               r = idx1(b);
            end
         else
            r = idx1;
         end
         % HKH, add break if r isempty
         if isempty(r), break; end

         S = [S setdiff(idx1,S)];
%xprinti(S,'S:',4,23);

         if f_star > F(r) % We must move horisontally to the left
            f_star = F(r);
         end

         % Compute the intersection of the line y(x) = slope*x + const,
         % with all curves not in [idx;r_previous]
         slope = -1/D(r);
         const = h(r)-slope*f_star;
 
         % Try speed up search for f_star_max just checking relevant rectangles

         %idxspeed = find( (cGD<cGD(r)) | (D>D(r)) );
         %idxspeed = setdiff(idxspeed,[idx1 r_previous]);  

         if m > 0
            ix0 =  (cGD<cGD(r)) | (D(1:nFunc)>D(r)) ;
         else
            ix0 =  D(1:nFunc)>D(r);
         end
         ix0(r_previous) = 0;
         % HKH: MUST SET THIS!
         ix0(ignoreIdx)  = 0; % not consider rects being fathomed
         ix0(idx1)       = 0;
         ix = find(ix0);

         %if isempty(ix) & isempty(idxspeed)
         %else
         %if ~all(ix == idxspeed)
         %   keyboard
         %end
         %end
         %f_star_m = f_star_max;

         if ~isempty(ix)
            if m > 0
               tmp1 = (cGD(ix)-const)/slope;
            else
               tmp1 = (-const/slope)*ones(1,length(ix));
            end

            iy = find(tmp1 < f_star & tmp1 >= F(ix));
            %f_star_max = max([f_star_max,tmp1(iy)]);
            tmpmax = max(tmp1(iy));
            if ~isempty(tmpmax)
               f_star_max = max(f_star_max,tmpmax);
            end
         end

         sl_hat = -1./D(ix);
         iz = find(sl_hat > slope);
         if m > 0
            con_hat = cGD(ix(iz)) - sl_hat(iz).*F(ix(iz));
         else
            con_hat = -sl_hat(iz).*F(ix(iz));
         end
         tmp2 = (con_hat-const)./(slope-sl_hat(iz));
         iv = find(tmp2 < f_star);
         if ~isempty(iv)
            % f_star_max = max([f_star_max,tmp2(iv)]);
            tmpmax = max(tmp2(iv));
            if ~isempty(tmpmax)
               f_star_max = max(f_star_max,tmpmax);
            end
         end

         %for dummy=1:length(idxspeed)
         %    i = idxspeed(dummy);
         %    % First assume that h(i) is constant i.e. max(F(i)-tmp1,0)=0
         %    tmp1 = (cGD(i)-const)/slope; % same as tmp1=-D(r)*(cGD(i)-const);
         %    if (tmp1 < f_star) & (tmp1 >= F(i) )
         %       % intersection between y(x) and h(i) where h(i) is constant.
         %       if tmp1 > f_star_max
         %          f_star_max = tmp1;
         %       end
         %    end
         %    %else % HKH-Maybe this should be an end statement !!!!!!!!
         %       % Assume h(i) is not constant i.e. max(F(i)-tmp2,0)=F(i)-tmp2
         %       slope_hat = -1/D(i);
         %       % Else, intersection not occurs or occur for values >= f_star
         %       if slope_hat > slope 
         %          %const_hat = (F(i)-tmp1)/D(i) + cGD(i) - slope_hat*f_star;
         %          const_hat = cGD(i) - slope_hat*F(i);
         %          tmp2 = (const_hat-const)/(slope-slope_hat);
         %          if tmp2 < f_star & tmp2 > f_star_max
         %             f_star_max = tmp2;
         %          end
         %       end
         %    %end
         %end % dummy=1:length(idxspeed)   

         %if f_star_max ~= f_star_m
         %   f_star_max
         %   f_star_m
         %   keyboard
         %end
 
         if isfinite(f_star_max) 
            f_star = f_star_max;
            r_previous = r;
         else
            % if curve r never intersected with another one
            break;
         end
         % DO NOT DOUBLE-COUNT ANY RECTANGLE !!!!!!   
      end % while 1  
   end % if ~feasible
   D(ixD) = 0; % Reset lengths to 0


   % -----------------------------------
   % STEP 3, CHOOSE ANY RECTANGLE IN S. 
   % -----------------------------------


%[Iter,length(S)]
%fprintf('Iter %d length(S) %d ',Iter,length(S))
%xprinti(S,'S:',4,23);

   for dummy=1:length(S) % for each rectangle in S
       r = S(dummy);


       % ---------------------------------------
       % STEP 4, TRISECT AND SAMPLE RECTANGLE r
       % ---------------------------------------

       if BRANCH == 0
          % CHOOSE A SPLITTING DIMENSION BY IDENTIFYING THE SET OF LONG SIDES OF
          % RECTANGLE r AND THEN CHOOSE THE LONG SIDE WITH THE SMALLEST t(i) VALUE.
          % IF MORE THAN ONE SIDE IS TIED FOR THE SMALLEST t(i) VALUE, CHOOSE THE 
          % ONE WITH THE LOWEST DIMENSIONAL INDEX.
          idx2 = find(Split(:,r)==min(Split(:,r)));
          if length(idx2) > 1
             [delta idx3] = min(t(idx2));
             i = idx2(idx3);
             % idx3 = find(t(idx2)==min(t(idx2)));
             %i = idx2(idx3(1));
             %%%i = idx2(idx3(end));
          else
             i = idx2;
          end
       elseif BRANCH == 1
          if ~isempty(IntVars)
             delta= min(Split(IntVars,r));
             if isinf(delta)
                i = [];
             else
                idx2 = find(Split(IntVars,r)==delta);
                if length(idx2) > 1
                   % Pick the first index with lowest t
                   [delta idx3] = min(t(IntVars(idx2)));
                   i = IntVars(idx2(idx3));
                else
                   i = IntVars(idx2);
                end
             end
          else
             i = [];
          end
          if isempty(i)
             idx2 = find(Split(:,r)==min(Split(:,r)));
             if length(idx2) > 1
                [delta idx3] = min(t(idx2));
                i = idx2(idx3);
             else
                i = idx2;
             end
          end
       elseif BRANCH == 2
          % First split rectangles with IntVars, range 1-3, i.e [0,1],[0,2] 
          if ~isempty(IntVars)
             delta= min(Split(IntVars(IX1),r));
             if isinf(delta)
                i = [];
             else
                idx2 = find(Split(IntVars(IX1),r)==delta);
                if length(idx2) > 1 % Pick the first index with lowest t
                   [delta idx3] = min(t(IntVars(IX1(idx2))));
                   i = IntVars(IX1(idx2(idx3)));
                else
                   i = IntVars(IX1(idx2));
                end
             end
          else
             i = [];
          end
          % Now split the other rectangles with IntVars
          if isempty(i)
             delta= min(Split(IntVars(IX2),r));
             if isinf(delta)
                i = [];
             else
                idx2 = find(Split(IntVars(IX2),r)==delta);
                if length(idx2) > 1
                   % Pick the first index with lowest t
                   [delta idx3] = min(t(IntVars(IX2(idx2))));
                   i = IntVars(IX2(idx2(idx3)));
                else
                   i = IntVars(IX2(idx2));
                end
             end
          end
          % Last split the continuous
          if isempty(i)
             idx2 = find(Split(:,r)==min(Split(:,r)));
             if length(idx2) > 1
                [delta idx3] = min(t(idx2));
                i = idx2(idx3);
             else
                i = idx2;
             end
          end
       elseif BRANCH == 3
          if ~isempty(IntVars)
             idx2 = find(~isinf(Split(IntVars,r)));
             if isempty(idx2)
                i = [];
             else
                delta = min(VarWeight(idx2));
                idx3 = find(VarWeight(idx2)==delta);
                if length(idx3) > 1
                   % Pick the first index with lowest t
                   [delta idx4] = min(t(IntVars(idx2(idx3))));
                   i = IntVars(idx2(idx3(idx4)));
                else
                   i = IntVars(idx2(idx3));
                end
             end
          else
             i = [];
          end
          if isempty(i)
             idx2 = find(Split(:,r)==min(Split(:,r)));
             if length(idx2) > 1
                [delta idx3] = min(t(idx2));
                i = idx2(idx3);
             else
                i = idx2;
             end
          end
       end

       % Updates
       t(i) = t(i)+1;

       %T(r) = T(r)+1;

       %% Update D(r)
       %j = mod(T(r),n-nFix);
       %k = (T(r)-j)/(n-nFix);
       %D(r) = 3^(-k)*sqrt(j/9+n-nFix-j)/2;

       %e_i    = [zeros(i-1,1);1;zeros(n-i,1)];
       e_i    = i;
       Split(i,r) = Split(i,r)+1;

       rightchild = 1; % flag if there will be a right/left child or not,
       % leftchild  = 1; % used when splitting along an integer dimension.

       % ******* LEFT NEW POINT ********
       % HKH. Move increment, do it directly first of all
       nFunc  = nFunc+1;
       Split_left = Split(:,r);
       if isMIP
          I_L_left = full(iL(:,r));
          I_U_left = full(iU(:,r));
       end   
       if I_logic(i) % We shall split along an integer dimension
          I_i   = find(IntVars==i);
          aa = iL(I_i,r);
          bb = iU(I_i,r);
          delta = floor( (bb-aa+1)/3 );
          c_left = full(C(:,r));
          if delta >= 1
             I_L_left(I_i) = aa;
             I_U_left(I_i) = aa+delta-1;
             iL(I_i,r) = aa + delta;
             iU(I_i,r) = bb - delta;
          elseif delta == 0 % Now there will be only 1 child. Left or right?
             %parent_center = (c_left(i)-x_L(i))./max(x_D(i),1E-20);
             parent_center = x_L(i) + c_left(i)*x_D(i);
             if abs(aa-parent_center) < 1E-4 % if aa==parent_center
                % leftchild  = 0;
                rightchild = 0;
                I_L_left(I_i) = bb;
                I_U_left(I_i) = bb;
                iL(I_i,r) = aa;
                iU(I_i,r) = aa;
             else
                rightchild=0;
                I_L_left(I_i) = aa;
                I_U_left(I_i) = aa;
                iL(I_i,r) = bb;
                iU(I_i,r) = bb;
             end
          else
             % rightchild = 0;
             disp('Error in glcSolve, this should not happen');
             return;
          end
          x_i_mid = floor((I_U_left(I_i)+I_L_left(I_i))/2);
          c_left(i) = (x_i_mid-x_L(i))/max(x_D(i),1E-20);
          if I_L_left(I_i)==I_U_left(I_i)
             Split_left(i) = Inf;
          end
          if iL(I_i,r)==iU(I_i,r)
             Split(i,r) = Inf;
          end
       else
          T(r) = T(r)+1;
          delta  = 3^(-Split(i,r));
          if delta < xTol
             Split(i,r)    = inf;
             Split_left(i) = inf;
          end
          %c_left = C(:,r) - delta*e_i; % Centerpoint for new left rectangle
          % Centerpoint for new left rectangle
          c_left = full(C(:,r));
          c_left(e_i) = c_left(e_i) - delta; 
       end
       if all(isinf(Split(:,r)))
          ignoreIdx = [ignoreIdx;r];
          cEND = 1;
       else
          cEND = 0;
       end
       if all(isinf(Split_left))
          ignoreIdx = [ignoreIdx;nFunc];
          lEND = 1;
       else
          lEND = 0;
       end
       % Update D(r)
       if nCont > 0
          j = mod(T(r),nCont);
          k = (T(r)-j)/(nCont);
       end
       if isMIP
          if nCont > 0
             z = (3^(-2*k))/4*(j/9+nCont-j);
          else
             z = 0;
          end
          % Length for mid point rectangle
          d = ceil(0.5*(iU(:,r) - iL(:,r)))./max(1,x_D(IntVars));
          D(r) = sqrt(z + sum(d.^2));
          %d = (iU(:,r) - iL(:,r));
          %d(d==1) = 2;    % Correct length 1 integer variables to length 1
          %D(r) = sqrt(z + sum(d.^2/4));

          % Length for left child rectangle
          d = ceil(0.5*(I_U_left - I_L_left))./max(1,x_D(IntVars));
          Dl = sqrt(z + sum(d.^2));
          %d = (I_U_left - I_L_left);
          %d(d==1) = 2;    % Correct length 1 integer variables to length 1
          %Dl = sqrt(z + sum(d.^2/4));

          %D  = [D Dl];
          D(nFunc)  = Dl;
          %iL = [iL I_L_left];
          %iU = [iU I_U_left];
          iL(:,nFunc) = I_L_left;
          iU(:,nFunc) = I_U_left;
       else
          D(r) = 3^(-k)*sqrt(j/9+nCont-j)/2;
          %D     = [D D(r)];
          D(nFunc)  = D(r);
       end
      
       % Transform c_left to original search space
       % x_left = x_L + c_left.*x_D;  
       x_left = tomsol(9,x_L, c_left,x_D); 
       if dLin > 0
          Ax_left = A*x_left;
          Afeas   = sum(max(0,max(b_L-bTol-Ax_left,Ax_left-bTol-b_U)));
       else
          Ax_left = [];
          Afeas   = 0;
       end   
       if dLin > 0 
          % Set bounds for center point rectangle
          cC      = full(C(:,r));
          cD      = 3.^(-Split(:,r));
          cL      = cC - 0.5*cD; 
          cU      = cC + 0.5*cD; 
          xL      = tomsol(9,x_L,cL,x_D); 
          xU      = tomsol(9,x_L,cU,x_D); 
          if isMIP
             xL(IntVars) = iL(:,r);
             xU(IntVars) = iU(:,r);
          end
          %xL    = max(x_L,xL);
          %xU    = min(x_U,xU);
       end
       if LinE > 0 
          % Test center point, if linear equality constraints feasible
          for k = 1:LinE
              j     = leq(k);
              Arow  = A(j,:);
              cKILL = CheckAxb(Arow,b_U(j),xL,xU);
              if cKILL
                 if ~cEND, ignoreIdx = [ignoreIdx;r]; end
                 break;
              end
          end
       else
          cKILL = 0;
       end
       % Avoid test if infeasible linear equality constraints
       if LinI > 0 & ~cKILL 
          % Test center point, if linear inequality constraints feasible 
          %inside rectangle
          for k = 1:LinI
              if k > length(b_L_idx)
                 j    = b_U_idx(k-length(b_L_idx));
                 Arow = A(j,:);
              else
                 j    = b_L_idx(k);
                 Arow = -A(j,:);
              end
              cKILL = CheckAxbU(Arow,g_U(k),xL,xU);
              if cKILL
                 if ~cEND, ignoreIdx = [ignoreIdx;r]; end
                 break;
              end
          end
       end
       if dLin > 0 
          % Set bounds for left point rectangle
          if I_logic(i) % We shall split along an integer dimension
             xL(e_i) = I_L_left(I_i);
             xU(e_i) = I_U_left(I_i);
          else
             Step = delta*x_D(e_i);
             % Just move bounds for selected rectangle side delta to the left
             xL(e_i) = xL(e_i) - Step;
             xU(e_i) = xU(e_i) - Step; 
          end
       end
       if LinE > 0 
          % Test left point, if linear equality constraints feasible
          for k = 1:LinE
              j     = leq(k);
              Arow  = A(j,:);
              lKILL = CheckAxb(Arow,b_U(j),xL,xU);
              if lKILL
                 if ~lEND, ignoreIdx = [ignoreIdx;nFunc]; end
                 break;
              end
          end
       else
          lKILL = 0;
       end

       if LinI > 0 & ~lKILL
          % Test left point, if non-feasible linear constraints 
          % are feasible inside rectangle

          gx_left = [-Ax_left(b_L_idx);Ax_left(b_U_idx)];
          tmpL    = gx_left - g_U(1:LinI);
          ix      =  find(tmpL(1:LinI) > g_T(1:LinI) ); % New point feasible ?

          for kk = 1:length(ix)
              k = ix(kk);
              if k > length(b_L_idx)
                 j    = b_U_idx(k-length(b_L_idx));
                 Arow = A(j,:);
              else
                 j    = b_L_idx(k);
                 Arow = -A(j,:);
              end
              lKILL = CheckAxbU(Arow,g_U(k),xL,xU);
              if lKILL
                 if ~lEND, ignoreIdx = [ignoreIdx;nFunc]; end
                 break;
              end
          end
       end
       if (fcALL > 0 | lKILL == 0) & ...
          (AxFeas == 0 | (AxFeas == 1 & (feasible|Afeas == 0|AvFunc == 0)))
          if SaveRes
             Prob.P = nFunc;
             [f_left,Res] = nlp_f(x_left, Prob, varargin{:});
          else
             f_left       = nlp_f(x_left, Prob, varargin{:});
          end
       else
          f_left = dBIG;
          dSum   = dSum +1;
       end
       if NLinI + NLinE > 0
          if fcALL > 0 | lKILL == 0 
             cx_left = nlp_c(x_left, Prob, varargin{:});
             cx_left = cx_left(:);
          else
             cx_left = dBIG*ones(NLinI+NLinE,1);
          end
       else
          cx_left = [];
       end
       if m > 0
          gx_left    = [-Ax_left(b_L_idx);Ax_left(b_U_idx); ...
                        -cx_left(c_L_idx);cx_left(c_U_idx); ... 
                         Ax_left(leq);cx_left(nleq)]; 
          tmpL = gx_left - g_U;
          if mE > 0
             rE = sum(beta.*abs(tmpL(mI+1:m)));
          else
             rE = 0;
          end

          G(:,nFunc) = tmpL; % Subtract g_U to get g(x)<=0
          % G(LinI+1:mI,nFunc) = tmpL(LinI+1:mI); % Subtract g_U to get g(x)<=0

          NowFeas = 0;
          if all( tmpL(1:mI) <= g_T(1:mI) ) % New point feasible ?
             if ~feasible
                feasible = 1;         % first feasible point
                if rE <= fEqualW 
                   glcfMin  = f_left; 
                else
                   glcfMin  = f_left + rE; 
                end
                NowFeas  = 1;
                fMinIdx  = nFunc;
                xBest    = x_left;
                SubRes   = Res;
                fMinEQ   = rE;
             elseif rE <= fEqualW & any(fMinEQ > fEqualW) 
                % Remove any nonfeasible equality solutions. 
                % Update glcfMin and fMinIdx
                ix      = find(fMinEQ <= fEqualW);
                if isempty(ix)
                   fBest   = Inf;
                else
                   [fBest,k] = min(F(fMinIdx(ix)));
                   fBest     = fBest + fMinEQ(k);
                end
                if isempty(ix) | f_left < fBest
                   glcfMin = f_left;
                   fMinIdx = nFunc;
                   xBest   = x_left;
                   SubRes  = Res;
                   fMinEQ  = rE;
                   glcfOld = Inf; % Must set Inf for IterPrint == 1
                elseif f_left+rE <= fBest
                   % Close point, add to set of minimum points
                   fMinIdx  = [nFunc,fMinIdx(ix)];
                   xBest    = [x_left,xBest(:,ix)];
                   SubRes   = Res;
                   fMinEQ   = [rE,fMinEQ(ix)];
                else
                   % Remove nonfeasible equality solutions
                   if PriLev > 0 | IterPrint > 0
                      'REMOVE left EQs'
                   end
                   fMinIdx  = fMinIdx(ix);
                   xBest    = xBest(:,ix);
                   fMinEQ   = fMinEQ(ix);
                end
             elseif rE > fEqualW & fMinEQ <= fEqualW  
                  % No test should be made, bad point
             elseif f_left + rE < glcfMin % Update glcfMin
                if glcfMin - (f_left + rE) < fEqualW
                   % Close point, add to set of minimum points
                   glcfMin = f_left + rE;
                   fMinIdx = [nFunc,fMinIdx];
                   xBest   = [x_left,xBest];
                   SubRes  = Res;
                   fMinEQ  = [rE,fMinEQ];
                else
                   if rE <= fEqualW 
                      glcfMin  = f_left; 
                   else
                      glcfMin  = f_left + rE; 
                   end
                   fMinIdx = nFunc;
                   xBest   = x_left;
                   SubRes  = Res;
                   fMinEQ  = rE;
                end
             elseif f_left + rE < glcfMin + fEqualW % Update fMinIdx
                % Close point, add to set of minimum points
                fMinIdx = [fMinIdx,nFunc];
                xBest   = [xBest,x_left];
                fMinEQ   = [fMinEQ,rE];
             end
          else
             InFeasL = find( tmpL(1:mI) > g_T(1:mI)); % Infeasible constraints
          end
       else
          if f_left < glcfMin % Update glcfMin
             if glcfMin - (f_left + rE) < fEqualW
                % Close point, add to set of minimum points
                fMinIdx = [nFunc,fMinIdx];
                xBest   = [x_left,xBest];
                SubRes  = Res;
                fMinEQ   = [rE,fMinEQ];
             else
                glcfMin = f_left;
                fMinIdx = nFunc;
                xBest   = x_left;
                SubRes  = Res;
                fMinEQ   = rE;
             end
          elseif f_left < glcfMin + fEqualW % Update fMinIdx
             % Close point, add to set of minimum points
             fMinIdx = [fMinIdx,nFunc];
             xBest   = [xBest,x_left];
             fMinEQ   = [fMinEQ,rE];
          end      
       end

       C(:,nFunc)     = c_left;
       F(nFunc)       = f_left;
       Split(:,nFunc) = Split_left;
       T(nFunc) = T(r);
      
       if ~feasible, fMin=Inf; else fMin=glcfMin; end
      
       InFeasR = [];

       if rightchild
          % ******* RIGHT NEW POINT ********
          % Move nFunc increment, do it first of all
          nFunc  = nFunc+1;
          Split_right = Split(:,r);
          if isMIP
             I_L_right = full(iL(:,r));
             I_U_right = full(iU(:,r));
          end   
          i = e_i; % Just for safety, reset i
          if I_logic(i) % We shall split along an integer dimension
             c_right = full(C(:,r));
             I_L_right(I_i) = bb-delta+1;
             I_U_right(I_i) = bb;
            
             x_i_mid = floor((I_U_right(I_i)+I_L_right(I_i))/2);
             %c_right(i) = (x_i_mid-x_L(i))/max(x_D(i),1E-20);
             % Not possible to get a rightchild if x_D(i) == 0
             c_right(i) = (x_i_mid-x_L(i))/x_D(i);
             %HKH-BUG in fortran???? if I_L_left(I_i)==I_U_left(I_i)
             if I_L_right(I_i)==I_U_right(I_i)
                Split_right(i) = Inf;
             end
         else
            %c_right = C(:,r) + delta*e_i; % Centerpoint for new right rectangle
            % Centerpoint for new right rectangle
            c_right      = full(C(:,r));
            c_right(e_i) = c_right(e_i) + delta;
         end
         if all(isinf(Split_right))
            ignoreIdx = [ignoreIdx;nFunc];
            rEND = 1;
         else
            rEND = 0;
         end
         % Transform c_right to original search space
         %x_right = x_L + c_right.*x_D;  
         x_right = tomsol(9,x_L, c_right,x_D); 
         if dLin > 0
            Ax_right = A*x_right;
            Afeas    = sum(max(0,max(b_L-bTol-Ax_right,Ax_right-bTol-b_U)));
            %format long
            %[b_L, Ax_right, b_U]
         else
            Ax_right = [];
            Afeas    = 0;
         end   
         if dLin > 0 
            % Set bounds for right point rectangle
            if I_logic(i) % Split along an integer dimension
               xL(e_i) = I_L_right(I_i);
               xU(e_i) = I_U_right(I_i);
            else
               Step = 2*delta*x_D(e_i);
               % Move bounds for selected rectangle side 2*delta to the right
               xL(e_i) = xL(e_i) + Step;
               xU(e_i) = xU(e_i) + Step; 
            end
         end
         if LinE > 0 
            % Test right point, if linear equality constraints feasible
            for k = 1:LinE
                j     = leq(k);
                Arow  = A(j,:);
                rKILL = CheckAxb(Arow,b_U(j),xL,xU);
                if rKILL
                   if ~rEND, ignoreIdx = [ignoreIdx;nFunc]; end
                   break;
                end
            end
         else
            rKILL = 0;
         end
         if LinI > 0 & rKILL == 0
            % Test right point, if non-feasible linear constraints 
            % are feasible inside rectangle

            gx_right = [-Ax_right(b_L_idx);Ax_right(b_U_idx)];
            tmpR     = gx_right - g_U(1:LinI);
            ix       = find(tmpR(1:LinI) > g_T(1:LinI) ); % New point feasible ?

            for kk = 1:length(ix)
                k = ix(kk);
                if k > length(b_L_idx)
                   j    = b_U_idx(k-length(b_L_idx));
                   Arow = A(j,:);
                else
                   j    = b_L_idx(k);
                   Arow = -A(j,:);
                end
                rKILL = CheckAxbU(Arow,g_U(k),xL,xU);
                if rKILL
                   if ~rEND, ignoreIdx = [ignoreIdx;nFunc]; end
                   break;
                end
            end
         end
         if (fcALL > 0 | rKILL == 0) & ...
            (AxFeas == 0 | (AxFeas == 1 & (feasible | Afeas == 0 | AvFunc == 0)))
            if SaveRes
               Prob.P = nFunc;
               [f_right,Res] = nlp_f(x_right, Prob, varargin{:});
            else
               f_right       = nlp_f(x_right, Prob, varargin{:}); 
            end
         else
            f_right = dBIG;
            dSum    = dSum +1;
         end
         if NLinI > 0 | ~isempty(nleq)
            if fcALL > 0 | rKILL == 0 
               cx_right = nlp_c(x_right, Prob, varargin{:});
               cx_right = cx_right(:);
            else
               cx_right = dBIG*ones(NLinI+NLinE,1);
            end
         else
            cx_right = [];
         end
         if m > 0
            gx_right    = [-Ax_right(b_L_idx);Ax_right(b_U_idx); ...
                           -cx_right(c_L_idx);cx_right(c_U_idx); ... 
                            Ax_right(leq);cx_right(nleq)]; 

            tmpR = gx_right - g_U;

            if mE > 0
               rE = sum(beta.*abs(tmpR(mI+1:m)));
            else
               rE = 0;
            end

            G(:,nFunc) = tmpR; % Subtract g_U to get g(x)<=0

            if all( tmpR(1:mI) <= g_T(1:mI) ) % New point feasible ?
               if ~feasible
                  feasible = 1;         % first feasible point
                  if rE <= fEqualW 
                     glcfMin  = f_right; 
                  else
                     glcfMin  = f_right + rE; 
                  end
                  NowFeas  = 1;
                  fMinIdx  = nFunc;
                  xBest    = x_right;
                  SubRes   = Res;
                  fMinEQ   = rE;
               elseif rE <= fEqualW & any(fMinEQ > fEqualW) 
                  % Remove any nonfeasible equality solutions. 
                  % Update glcfMin and fMinIdx
                  ix      = find(fMinEQ <= fEqualW);
                  if isempty(ix)
                     fBest   = Inf;
                  else
                     [fBest,k] = min(F(fMinIdx(ix)));
                     fBest     = fBest + fMinEQ(k);
                  end
                  if isempty(ix) | f_right < fBest
                     glcfMin = f_right;
                     fMinIdx = nFunc;
                     xBest   = x_right;
                     SubRes  = Res;
                     fMinEQ  = rE;
                     glcfOld = Inf; % Must set Inf for IterPrint == 1
                  elseif f_right+rE <= fBest
                     % Close point, add to set of minimum points
                     fMinIdx  = [nFunc,fMinIdx(ix)];
                     xBest    = [x_right,xBest(:,ix)];
                     SubRes   = Res;
                     fMinEQ   = [rE,fMinEQ(ix)];
                  else
                     % Remove nonfeasible equality solutions
                     if PriLev > 0 | IterPrint > 0
                        'REMOVE right EQs'
                     end
                     fMinIdx  = fMinIdx(ix);
                     xBest    = xBest(:,ix);
                     fMinEQ   = fMinEQ(ix);
                  end
               elseif rE > fEqualW & fMinEQ <= fEqualW  
                  % No test should be made, bad point
               elseif f_right + rE < glcfMin % Update glcfMin
                  if glcfMin - (f_right + rE) < fEqualW
                     % Close point, add to set of minimum points
                     glcfMin  = f_right + rE;
                     fMinIdx  = [nFunc,fMinIdx];
                     xBest    = [x_right,xBest];
                     SubRes   = Res;
                     fMinEQ   = [rE,fMinEQ];
                  else
                     if rE <= fEqualW 
                        glcfMin  = f_right; 
                     else
                        glcfMin  = f_right + rE; 
                     end
                     fMinIdx  = nFunc;
                     xBest    = x_right;
                     SubRes   = Res;
                     fMinEQ   = rE;
                  end
               elseif f_right + rE < glcfMin + fEqualW % Update fMinIdx
                  % Close point, add to set of minimum points
                  fMinIdx = [fMinIdx,nFunc];
                  xBest   = [xBest,x_right];
                  fMinEQ   = [fMinEQ,rE];
               end
            else
              InFeasR = find( tmpR(1:mI) > g_T(1:mI)); % Infeasible constraints
            end
         elseif f_right < glcfMin % Update glcfMin and x_min
            if glcfMin - (f_right + rE) < fEqualW
               % Close point, add to set of minimum points
               fMinIdx  = [nFunc,fMinIdx];
               xBest    = [x_right,xBest];
               SubRes   = Res;
               fMinEQ   = [rE,fMinEQ];
            else
               glcfMin  = f_right;
               fMinIdx  = nFunc;
               xBest    = x_right;
               SubRes   = Res;
               fMinEQ   = rE;
            end
         elseif f_right < glcfMin + fEqualW % Update fMinIdx
            % Close point, add to set of minimum points
            fMinIdx = [fMinIdx,nFunc];
            xBest   = [xBest,x_right];
            fMinEQ   = [fMinEQ,rE];
         end   
         if isMIP
            %iL = [iL I_L_right];
            %iU = [iU I_U_right];
            iL(:,nFunc) = I_L_right;
            iU(:,nFunc) = I_U_right;
            % Length for right child rectangle
            d = ceil(0.5*(I_U_right - I_L_right))./max(1,x_D(IntVars));
            Dr = sqrt(z + sum(d.^2));
            %d = (I_U_right - I_L_right);
            %d(d==1) = 2;    % Correct length 1 integer variables to length 1
            %Dr = sqrt(z + sum(d.^2/4));
            %D    = [D Dr];
            D(nFunc)  = Dr;
         else
            %D     = [D D(r)];
            D(nFunc)  = D(r);
         end
         C(:,nFunc)     = c_right;
         F(nFunc)       = f_right;
         Split(:,nFunc) = Split_right;
         T(nFunc) = T(r);
         
         if ~feasible, fMin=Inf; else fMin=glcfMin; end

      end % if rightchild

      if PriLev > 0 | IterPrint > 0
         if NowFeas > 0
            fprintf('Feasible. Center %d Left/Right %d %f\n',r,nFunc,F(r));
         end
      end
      if AxFeas == 1 & NowFeas == 1
         % Recompute f(x) for center, if not computed
         if F(r) >= dBIG
            % Use x_right as temp storage
            x_right  = tomsol(9,x_L,cC,x_D); 
            if SaveRes
               Prob.P = r;
               [F(r),Res] = nlp_f(x_right, Prob, varargin{:});
            else
               F(r)       = nlp_f(x_right, Prob, varargin{:}); 
            end
            % fprintf('Recompute %d %f\n',r,F(r));
         end
      end
   
      if m > 0 & useRoC > 0
         % UPDATE THE s(j):s
         % Transform C(:,r) to original search space
         %x_mid = x_L + C(:,r).*x_D;  
         x_mid = tomsol(9,x_L, full(C(:,r)),x_D); 
         norm_left  = norm(x_left-x_mid);
         if AvFunc < AvIter
            if norm_left ~= 0 & f_left<dBIG & (fcALL==2 | (lKILL==0 & F(r)<dBIG))
               s_0 = (AvFunc*s_0 + abs(f_left-F(r))/norm_left)/(AvFunc+1);
               if NLinI > 0
                  ll = LinI+1;
                  s(ll:mI) = (AvFunc*s(ll:mI) + abs(tmpL(ll:mI) ...
                              - G(ll:mI,r))/norm_left)/(AvFunc+1);
               end
               if NLinE > 0
                  ll = mI+LinE+1;
                  s(ll:m)  = (AvFunc*s(ll:m) + abs( tmpL(ll:m) ...
                              - G(ll:m,r)) /norm_left)/(AvFunc+1);
               end
               AvFunc = AvFunc + 1;
            end
            if rightchild & f_right<dBIG & (fcALL==2 | (rKILL==0 & F(r)<dBIG))
               norm_right = norm(x_right-x_mid);
               if norm_right ~= 0 
                  s_0 = (AvFunc*s_0 + abs(f_right-F(r))/norm_right)/(AvFunc+1);
                  if NLinI > 0
                     ll = LinI+1;
                     s(ll:mI) = (AvFunc*s(ll:mI)+abs(tmpR(ll:mI) ...
                                 - G(ll:mI,r))/norm_right)/(AvFunc+1);
                  end
                  if NLinE > 0
                     ll = mI+LinE+1;
                     s(ll:m)  = (AvFunc*s(ll:m)+abs(tmpR(ll:m) ...
                                 - G(ll:m,r))/norm_right)/(AvFunc+1);
                  end
                  AvFunc = AvFunc + 1;
               end
            end
         else
            % Use exponential smoothing
            if norm_left~=0 & f_left<dBIG & (fcALL==2|(lKILL == 0 & F(r)<dBIG))
               s_0         = alfa*s_0 + (1-alfa)*abs(f_left-F(r))/norm_left;
               if NLinI > 0
                  ll = LinI+1;
                  s(ll:mI) = alfa*s(ll:mI) + (1-alfa) * ...
                             abs(tmpL(ll:mI)-G(ll:mI,r))/norm_left;
               end
               if NLinE > 0
                  ll = mI+LinE+1;
                  s(ll:m)  = alfa*s(ll:m) + (1-alfa) * ...
                             abs(tmpL(ll:m)-G(ll:m,r))/norm_left;
               end
               AvFunc = AvFunc + 1;
            end
            if rightchild & f_right<dBIG & (fcALL==2|(rKILL == 0 & F(r)<dBIG))
               norm_right  = norm(x_right-x_mid);
               if norm_right ~= 0 
                  s_0        = alfa*s_0 + (1-alfa)*abs(f_right-F(r))/norm_right;
                  if NLinI > 0
                     ll = LinI+1;
                     s(ll:mI)= alfa*s(ll:mI) + (1-alfa) * ...
                               abs(tmpR(ll:mI)-G(ll:mI,r))/norm_right;
                  end
                  if NLinE > 0
                     ll = mI+LinE+1;
                     s(ll:m) = alfa*s(ll:m) + (1-alfa) * ...
                               abs(tmpR(ll:m)-G(ll:m,r))/norm_right;
                  end
                  AvFunc = AvFunc + 1;
               end
            end
         end
         if PriLev > 2 
            fprintf('AvF %d s_0 %f',AvFunc,s_0)
            if LinI > 0
               fprintf(' Lin ')
               fprintf(' %f',s(1:LinI))
            end
            if NLinI > 0
               fprintf(' NonLin ')
               fprintf(' %f',s(LinI+1:mI))
            end
            if LinE > 0
               fprintf(' LinE')
               fprintf(' %f',s(mI+1:mI+LinE))
            end
            if NLinE > 0
               fprintf(' NonLinE')
               fprintf(' %f',s(mI+LinE+1:m))
            end
            fprintf('\n')
         end
      end
      if PriLev > 1
         fprintf('It:%5d ',Iter);
         if nFuncOld > 0
            fprintf('Ev:%5d TotRec:%7d ',nFunc-nFuncOld,nFunc-nFuncOld+nRectTot);
         else
            fprintf('Ev:%5d ',nFunc);
         end
         fprintf('f:%16.10f',fMin);
         if rightchild
            fprintf(' %16.10f %16.10f',F(nFunc-1:nFunc));
         else
            fprintf(' %16.10f                 ',F(nFunc));
         end
         fprintf(' Var %2d',e_i);
         fprintf(' D:%8.4g',delta);
         if ~(dLin == 0 & isempty(IntVars))
            fprintf(' Ig:%4d',length(ignoreIdx));
         end
         fprintf(' Split %4d',Split(e_i,r));
         if rightchild
            fprintf(' %4d %4d',Split(e_i,nFunc-1:nFunc));
         else
            fprintf(' %4d    ',Split(e_i,nFunc));
         end
         if ~isempty(IntVars)
            fprintf(' IV:')
            for i = 1:length(IntVars)
                fprintf(' %d',x_left(IntVars(i))); 
            end
            if rightchild
               fprintf(' IV-r:')
               for i = 1:length(IntVars)
                   fprintf(' %d',x_right(IntVars(i))); 
               end
            end
         end
         if isinf(fMin)
            fprintf(' Feas %d',feasible);
            if m > 0, fprintf(' c:%20.10f',min(max(G(:,1:nFunc)))); end
         end
         fprintf('\n')
      end
      
   end % for each rectangle in S
   
   %
   %  STEP 5, UPDATE S (this step is handled by the for loop)
   %
   
   %
   %  STEP 6, ITERATE (this step is handled by the for loop)
   %
   
   if isempty(S) % If no feasible integer solution exist
      break;
   else
      
   if PriLev > 0 | (IterPrint > 0 & glcfMin < glcfOld)   
      fprintf('It:%5d ',Iter);
      if nFuncOld > 0
         fprintf('Ev:%5d TotRec:%7d ',nFunc-nFuncOld,nFunc-nFuncOld+nRectTot);
      else
         fprintf('Ev:%5d ',nFunc);
      end
      fprintf('f:%16.10f',fMin);
      %fprintf(' D:%8.4g',delta);
      if ~(dLin == 0 & isempty(IntVars))
         fprintf(' Ig:%4d',length(ignoreIdx));
      end
      %fprintf(' Var %2d',e_i);
      %fprintf(' Split %3d',Split(e_i,nFunc));
      if mI > 0
         fprintf(' IE:%9.6f',max(full(G(:,fMinIdx(end)))));
      end
      if mE > 0
         fprintf(' Eq:%9.6f',fMinEQ(1));
      end
      fprintf(' e_i:%d',e_i);
      if isinf(fMin)
         if rightchild
            if ~isempty(IntVars) 
               fprintf(' IV:')
               fprintf(' %d',full(C(IntVars,nFunc-1))); 
            end
            fprintf(' ');
            fprintf(' %f',full(C(find(~I_logic),nFunc-1))); 
         end
         if ~isempty(IntVars) 
            fprintf(' IV:')
            fprintf(' %d',full(C(IntVars,nFunc))); 
         end
         fprintf(' ');
         fprintf(' %f',full(C(find(~I_logic),nFunc))); 
      else
         if ~isempty(IntVars) 
            fprintf(' IV:')
            for i = 1:length(IntVars)
                fprintf(' %d',xBest(IntVars(i),1)); 
            end
         end
         fprintf(' ');
         fprintf(' %f',xBest(find(~I_logic),1)); 
      end
      if isinf(fMin)
         %fprintf(' Feas %d',feasible);
         fprintf(' InFeasL');
         fprintf(' %d',InFeasL);
         if ~isempty(InFeasR)
            fprintf(' InFeasR');
            fprintf(' %d',InFeasR);
         end
         if m > 0, fprintf(' c:%13.10f',min(max(full(G(:,1:nFunc))))); end
      end
      fprintf('\n');
      glcfOld = glcfMin;
   end

   if convFlag == 0 & feasible
      convFlag = isClose(fGoal,fMin,fTol,nFunc,nFuncOld,nFuncTot,Iter,EpsGlob,PriLev);
   end
   
   end
   
end % while 1  (main loop)
   
%fprintf('\n\n');

% SAVE RESULTS
C     = C(:,1:nFunc);
Split = Split(:,1:nFunc);
F     = F(1:nFunc);
T     = T(1:nFunc);
D     = D(1:nFunc);
if isMIP
   iL = iL(:,1:nFunc);
   iU = iU(:,1:nFunc);
end
if m > 0
   % Save C in transformed form
   G  = G(:,1:nFunc);
end
Result = ResultDef(Prob);
Result.Solver = 'glcSolve';
Result.SolverAlgorithm = 'Constrained DIRECT (Jones 2001)';

% Maximal size of possible triangles
tmpset1=ones(1,nFunc);
tmpset1(ignoreIdx) = 0; % not consider rects being fathomed

Result.maxTri = 3.^-(min(min(Split(:,find(tmpset1))))); 

if feasible
   Result.f_k  = glcfMin;    % Best function value
else
   % Result.f_k  = Inf;      % No feasible point found
   Result.f_k  = glcfMin;    % Best function value
end
Result.Iter = Iter;     % Number of iterations

% For restart

% Find all points i with F(i)=glcfMin
if feasible
   Result.minPnts = fMinIdx;
   Result.x_k     = xBest;
   if convFlag == 0
      if Iter >= MaxIter 
         Result.Inform  = 3;
      elseif nFunc >= MaxFunc
         Result.Inform  = 4;
      else
         Result.Inform  = 0;
      end
   else
      Result.Inform  = convFlag;
   end

   %idx = find(F==glcfMin);
   %if m > 0 % If there are constraints, pick out the feasible points in idx.
   %   if mE == 0
   %      idx2 = all(G(:,idx) < cTol ); % if feasible
   %   elseif mI == 0
   %      idx2 = all(abs(G(:,idx)) < cTol ); % if feasible
   %   else
   %      idx2 = all(G(1:mI,idx) < cTol) & all(abs(G(mI+1:m,idx)) < cTol);
   %   end
   %   idx = idx(idx2);
   %end
   %% All points i with F(i)=glcfMin, and feasible
   %Result.x_k = tomsol(9,x_L,C(:,idx),x_D);    
else
   Result.Inform  = 91;
   %LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
   %NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
   %mI    = LinI + NLinI;                     % Number of constraints
   %LinE  = length(leq);
   %NLinE = length(nleq);
   %mE    = LinE + NLinE;
   %m     = mI + mE;
   %g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
   %        b_L(leq);c_L(nleq)];
   %g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
   %       cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
   %       bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
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
      Fpen = [G(1:mI,:) - g_T(1:mI)*ones(1,N); ...
              abs(G(mI+1:m,:))-g_T(mI+1:m)*ones(1,N)];
   else
      Fpen = max([G(1:mI,:) - g_T(1:mI)*ones(1,N); ...
              abs(G(mI+1:m,:))-g_T(mI+1:m)*ones(1,N)]);
   end
   if Afeas == 1
      ixx = find(ix);
      [fBest,iMin]=min(F(ixx)+Fpen(ixx));
      fMinIdx = ixx(iMin);
   elseif Afeas == 2
      ixx = find(ix);
      [fBest,iMin]=min(Fpen(ixx));
      fMinIdx = ixx(iMin);
   elseif Afeas == 3
      [fBest,fMinIdx]=min(F+Fpen);
   else
      [fBest,fMinIdx]=min(Fpen);
   end
   Result.minPnts = fMinIdx;
   xBest          = tomsol(9,x_L,full(C(:,fMinIdx)),x_D);    
   Result.x_k     = xBest;
   Result.f_k     = F(fMinIdx);
   if Result.f_k == dBIG
      if SaveRes
         Prob.P = nFunc + 1;
         [Result.f_k,SubRes] = nlp_f(xBest, Prob, varargin{:});
      else
         Result.f_k          = nlp_f(xBest, Prob, varargin{:});
      end
   end
end

Result.SubResult = SubRes;
Result.FuncEv    = nFunc-dSum-nFuncOld;
nFuncTot         = nFunc-dSum-nFuncOld+nFuncTot;
nRectTot         = nFunc-nFuncOld+nRectTot;
Result.FuncEvTot = nFuncTot;

if NLinI + NLinE > 0
   c_k=[];
   for i=1:length(Result.minPnts)
      cxx=nlp_c(Result.x_k(:,i), Prob, varargin{:});
      c_k = [c_k  cxx(:)];
   end
   Result.c_k = c_k; % Constraint value at x_k;
   Result.ConstrEv = Result.FuncEv;
end
if LinI + LinE > 0 & size(Result.x_k,2) == 1
   Result.Ax=A*Result.x_k; % Linear constraint value at all x_k;
end
m                = nFunc;
if feasible
   Result.ExitFlag = 0;
   if Result.Inform==3
      Result.ExitText = ['Maximum Iterations reached. ' ...
                       'Tried ' num2str(nRectTot) ' rectangles, computing '...
                       num2str(nFuncTot) ' f(x) values'];
   else
      Result.ExitText = ['Tried ' num2str(nRectTot) ' rectangles, computing '...
                       num2str(nFuncTot) ' f(x) values'];
   end
else
   Result.ExitFlag = 7;
   Result.ExitText = ['Tried ' num2str(nRectTot) ' rectangles, computing '...
                       num2str(nFuncTot) ' f(x) values, no feasible point'];
end
if cpumax
   Result.Inform=9;
   Result.ExitText= [Result.ExitText '. Max CPU reached. '];
end

if ( PriLev >= 1 | IterPrint > 0 ) & dSum > 0 
   fprintf('Infeasible linear constraints for ');
   fprintf('%d out of %d calls, net calls %d\n',dSum,nFunc-nFuncOld,nFunc-nFuncOld-dSum);
end

WarmStartInfo.Name      = Name; 
WarmStartInfo.C         = C;  
WarmStartInfo.F         = F;  
WarmStartInfo.T         = T;  
WarmStartInfo.D         = D;  
WarmStartInfo.G         = G;  
WarmStartInfo.iL        = iL;  
WarmStartInfo.iU        = iU;  
WarmStartInfo.s_0       = s_0;  
WarmStartInfo.s         = s;  
WarmStartInfo.t         = t;  
WarmStartInfo.ignoreIdx = ignoreIdx;  
WarmStartInfo.feasible  = feasible; 
WarmStartInfo.Split     = Split;  
WarmStartInfo.glcfMin   = glcfMin; 
WarmStartInfo.fMinIdx   = fMinIdx; 
WarmStartInfo.fMinEQ    = fMinEQ; 
WarmStartInfo.SubRes    = SubRes; 
WarmStartInfo.nFuncTot  = nFuncTot; 
WarmStartInfo.nRectTot  = nRectTot;   

SaveWarmStartFile(WarmStartInfo);
Result.glcSolve.WarmStartInfo = WarmStartInfo;

Result=endSolve(Prob,Result);

function convFlag = isClose(fGoal,f,fTol,nFunc,nFuncOld,nFuncTot,Iter,EpsGlob,PriLev)

convFlag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convFlag = 1;
elseif fGoal == 0
   if abs(f-fGoal) < fTol
      convFlag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convFlag = 3;
end

if convFlag > 0 & PriLev >= 0 
   if convFlag == 1
      fprintf('\n\nFunction value %f is less than fGoal %f \n',f,fGoal);
   elseif convFlag == 2
      fprintf('\n\nError in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal),fTol);
   elseif convFlag == 3
      fprintf('\n\nRelative error in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('Number of rectangles:  %d\n',nFunc-nFuncOld);
   if nFuncOld > 0
      fprintf(' (In total %d)',nFunc-nFuncOld+nFuncTot);
   end
   fprintf('\n');
   fprintf('Number of iterations:            %d\n',Iter);
   fprintf('Epsilon:                         %f\n',EpsGlob);
end
if convFlag == 3, convFlag = 2; end

% ------------------------------------------------------------------------------------
function KILL = CheckAxbU(A,b,x_L,x_U)
% ------------------------------------------------------------------------------------
% A  Row in A matrix
% b  Right hand side value for the corresponding A row
% x_L Lower bounds on x in box
% x_U Upper bounds on x in box
%
% CheckAxbU returns KILL = 1 if Ax <= b is not possible to fulfill in the box

KILL  = 0;
x     = zeros(length(A),1);
ix    = find(A > 0);
x(ix) = x_L(ix);
ix    = find(A < 0);
x(ix) = x_U(ix);
if (A*x-b) > 1E-10
   % Not possible to get feasible inside rectangle
   % fprintf('CheckAxbU: Constraint not feasible in rectangle\n');
   KILL = 1;
end
%if all([1 0 0]'==x(1:3))
%   [x_L, x_U]
%   KILL
%end

% ------------------------------------------------------------------------------------
function KILL = CheckAxb(A,b,x_L,x_U)
% ------------------------------------------------------------------------------------
% A   Row in A matrix
% b   Right hand side value for the corresponding A row
% x_L Lower bounds on x in box
% x_U Upper bounds on x in box
%
% CheckAxb returns KILL = 1 if Ax == b is not possible to fulfill in the box

KILL  = 0;
x     = zeros(length(A),1);
iP    = find(A > 0);
x(iP) = x_L(iP);
iN    = find(A < 0);
x(iN) = x_U(iN);
Ax    = A*x;
Ax1 = Ax;
x1 = x;
if abs(Ax - b) <= 1E-10, return, end
POS   = Ax > b;
x(iP) = x_U(iP);
x(iN) = x_L(iN);
Ax    = A*x;
if (POS & (Ax-b) > 1E-10) | (~POS & (b-Ax) > 1E-10)
   % Not possible to get feasible inside rectangle
   %fprintf('CheckAxb: Constraint not feasible in rectangle\n');
   KILL = 1;
else
   %[x_L, x_U, x]
end
%KILL = 0;
%if KILL
%   A
%   [x_L x_U x1 x]
%   [Ax Ax1 b]
%end
%if all([1 0 0]'==x(1:3))
%   [x_L, x_U, x]
%   KILL
%keyboard
%end

function SaveWarmStartFile(WarmStartInfo)

Name      = WarmStartInfo.Name;
C         = WarmStartInfo.C;
F         = WarmStartInfo.F;
T         = WarmStartInfo.T;
D         = WarmStartInfo.D;
G         = WarmStartInfo.G;
iL        = WarmStartInfo.iL;
iU        = WarmStartInfo.iU;
s_0       = WarmStartInfo.s_0;
s         = WarmStartInfo.s;
t         = WarmStartInfo.t;
ignoreIdx = WarmStartInfo.ignoreIdx;
feasible  = WarmStartInfo.feasible;
Split     = WarmStartInfo.Split;
glcfMin   = WarmStartInfo.glcfMin;
fMinIdx   = WarmStartInfo.fMinIdx;
fMinEQ    = WarmStartInfo.fMinEQ;
SubRes    = WarmStartInfo.SubRes;
nFuncTot  = WarmStartInfo.nFuncTot;
nRectTot  = WarmStartInfo.nRectTot; 

try
   save('glcSave.mat','Name','C','F','T','D','G','iL','iU','s_0','s','t',...
      'ignoreIdx','feasible','Split','glcfMin','fMinIdx','fMinEQ',...
      'SubRes','nFuncTot','nRectTot');
catch
   warning('Failed to save warm start information to glcSave.mat');
   disp(lasterr);
   disp('Warm start information is available in Result.glcSolve.WarmStartInfo');
end

% MODIFICATION LOG:
%
% 990408  mbk  Modified to make restart possible.
% 990416  mbk  Output if PriLev > 0.
% 990416  mbk  Small changes in comments.
% 000830  hkh  Some speedups
% 000928  hkh  Revision for v3.0
% 001011  hkh  Include Name in glcSave.mat for safety
% 010330  hkh  Check for if ~isempty(tmpidx), row 449, to avoid 6.0 warning
% 010416  hkh  Define isClose, and target value convergence, tomsol speedups 
% 010416  hkh  Bug in speedup of the idxspeed computation
% 011031  hkh  Improve comments. Add maxTri calculation.
% 011110  hkh  Fixed errors in comments
% 011212  hkh  I(tmpidx) should be IUL(tmpidx), in error printout
% 020110  hkh  Change name f_min to glcfMin, f_minn to fMin
% 020309  hkh  Clean up code. New way of computing s_0, rateOfChange
% 020311  hkh  Fixed old bug computing rateOfChange. Weight linear/nonlinear
% 020313  hkh  Major revision, speedup of bottle necks.
% 020313  hkh  Change logic to check both slope and constant function
% 020325  hkh  Major revision
% 020325  hkh  Add LinWeight (linear /nonlinear balance) parameter
% 020325  hkh  Add fEqual parameter, tolerance when points are equal
% 020420  hkh  Incorrect check, checking left. not right int triangle
% 020506  hkh  Improving comments about DIRECT algorithm
% 030124  hkh  Allow general IntVars input similar to mipSolve
% 030124  hkh  Do preSolve if any lower or upper bound is inf or empty
% 030124  hkh  Accept empty lower or upper bound as input
% 031123  hkh  Add nFuncOld to fMinIdx on many places
% 031123  hkh  Change declaration of T if WarmStart. Other minor changes
% 040111  hkh  Change call to inisolve
% 040308  hkh  Change ignoreidx to ignoreIdx,I_L to iL,I_U to iU, IL to IUL
% 040308  hkh  Make ignoreidx column vector
% 040308  hkh  Avoiding calls if linear constraints infeasible
% 040328  hkh  Only check convergence (if close to fGoal) if feasible
% 040328  hkh  Add Inform value, revised ExitFlag values
% 040331  hkh  If infeasible, find point close to feasible
% 040415  hkh  Change beta to always be positive, doing abs()
% 040418  hkh  Save SubRes, additional output from suboptimization
% 040418  hkh  New input Prob.SaveRes, if 1 nlp_f returns 2nd output SubRes
% 040418  hkh  New output field Result.SubResult
% 040421  hkh  Declare arrays, to make faster, and run more iterations 
% 040422  hkh  Sparse or dense version, dependent on Prob.LargeScale
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 041017  hkh  Add input xIP and fIP, Use fIP when computing f_star
% 041123  hkh  Change call to tomRun
% 050402  hkh  Move nFunc increment to be first, always use nFunc as current
% 050402  hkh  Move D(r) computation before left point handling
% 050402  hkh  New functions CheckAxbU and CheckAxb for linear constraints
% 050405  hkh  Fathomed rectangles could pop up again, ix0(ignoreIdx)=0 added
% 050405  hkh  Move LinWeight, fEqual from Prob.GO to Prob.glcDirect
% 050405  hkh  Added iteration print of IE, max of inequality constraints
% 050405  hkh  Major revision. Add fcALL, AxFeas options
% 050405  hkh  New input useRoC, to use or not use rate of change
% 050405  hkh  New EqConFac, weight factor for equalities, 0 means all ones
% 050406  hkh  Squeeze arrays if WarmStart, delete all ignoreIdx rectangles
% 050406  hkh  Added variables nRectTot, nFuncTot in mat file, changed ExitText
% 050406  hkh  Made distinction clearer between rectangles and function values
% 050406  hkh  Add different branch strategies for MIP, input BRANCH
% 050406  hkh  Add different strategies RECTIE to handle ties, input RECTIE 
% 050407  hkh  Use relative tolerance fEqualW = fEqual*sum(beta)
% 050407  hkh  FSTAR, three different options how to compute target f_star
% 050412  hkh  Correction of xL,xU rectangle computation, also more efficient
% 050412  hkh  Must not delete rectangles, if in fMinIdx
% 050412  hkh  Add one new type of RECTIE, use distance D in two of them
% 050412  hkh  G was not squeezed after WarmStart, when nonempty(ignoreIdx)
% 050412  hkh  If LargeScale [], set 1 if MaxFunc >= 10000
% 050419  hkh  xL,xU rectangle incorrect for integers, use I_i and I_logic
% 050419  hkh  Set RECTIE=1 as default, still RECTIE=0 is sometimes faster
% 050419  hkh  Additional abs in some RoC computations, removed
% 050419  hkh  Add input forgetting factor alpha, default 0.9, alfa in code
% 050419  hkh  Add input AvIter, how long filter startup before exp.smoothing
% 050419  hkh  Test for rE<=fEqualW if 1st time feasible, when setting glcfMin
% 050420  hkh  Added new test to avoid f+rE being set if rE <=fEqualW
% 051006  hkh  Safe handling of MaxFunc, new default
% 060327  hkh  Add output of Result.ConstrEv=Result.FuncEv, if NLinI+NLinE > 0
% 060814  med  FUNCS used for callbacks instead
% 061003  ango Safe try-catch saving to glcSave.mat
% 061117  hkh  Improve ExitText if feasible
% 070222  hkh  Revise IntVars handling, use new format
