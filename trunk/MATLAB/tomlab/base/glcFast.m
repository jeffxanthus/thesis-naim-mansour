% glcFast.m
%
% glcSolve solves general constrained mixed integer global optimization
% This version, glcFast, is a Fortran MEX version of glcSolve
%
% glcSolve / glcFast implements the DIRECT algorithm by Donald R. Jones
% presented in the paper "DIRECT", Encyclopedia of Optimization,
% Kluwer Academic Publishers, 2001.
% The algorithm is expanded to handle nonlinear and linear equalities, and
% linear inequalities.
%
% glcFast solves problems of the form:
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
% function Result = glcFast(Prob)
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
%   PriLevOpt Print Level
%             0 = silent.
%             1 = Warm Start info, and same output as IterPrint = 1 (see below)
%             2 = Printing each iteration, Print iteration number, number of
%                 evaluated points, best f(x) so far and its x-values
%                 each iteration, where the best point has improved
%             3 = As PriLevOpt=2, also print sampled x, f(x) and if feasible
%
%   WarmStart If true, >0, glcFast reads the output from the last run from the
%             mat-file glcFastSave.mat, and continues from the last run.
%
%   MaxCPU    Maximum CPU Time (in seconds) to be used.
%
%
% ---------------------------------------
% glcDirect   Structure with DIRECT algorithm specific parameters. Fields used:
% ---------------------------------------
%
%
%   fcAll     =0 (Default). If linear constraints cannot be feasible anywhere
%             inside rectangle, skip f(x) and c(x) computation for middle point
%             =1. Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Do not update rates of change 
%             for the constraints.
%             =2. Always compute f(x) and c(x), even if linear constraints are
%             not feasible anywhere in rectangle. Update rates of change
%             constraints.
%
%   AxFeas    Set nonzero to make glcSolve skip f(x) evaluations, when the
%             linear constraints are infeasible, and still no feasible point
%             has been found.  The default is 0. Value 1 demands fcALL == 0.
%             This option could save some time if f(x) is a bit costly, however
%             overall performance could on some problems be dramatically worse
%           
%   fEqual    All points with function values within tolerance fEqual are
%             considered to be global minima and returned
%   LinWeight RateOfChange = LinWeight*|a(i,:)| for linear constraints.
%             Balance between linear and nonlinear constraints.  Default 0.1. 
%             The higher value, the less influence from linear constraints
%
% ---------------------------------------
% optParam    Structure in Prob, Prob.optParam.
% ---------------------------------------
%             Defines optimization parameters. Fields used:
%   MaxIter   Maximal number of iterations, default max(5000,n*1000);
%   MaxFunc   Maximal number of function evaluations, default max(10000,n*2000)
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
%             It is advised to number the integer values as the first
%             variables, before the continuous. The tree search will then
%             be done more efficiently.
%
%   fIP       An upper bound on the optimal f(x) value.  If empty, set as Inf.
%   xIP       The x-values giving the fIP value.
%             If fIP empty and xIP given, fIP will be computed
%             if xIP nonempty, its feasibility is checked
%
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  c_k      Nonlinear constraints values at x_k
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  maxTri   Maximum size of any triangle
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag 0 = Normal termination, max number of iterations /func.evals reached
%           2 = Some upper bounds below lower bounds
%           4 = Numerical trouble, and cannot continue
%           7 = Reached maxFunc or maxIter, NOT feasible
%           8 = Empty domain for integer variables
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
%               2*sLen .GE. maxFDim - nFunc, save warm start info, stop
%           7 = Space is dense
%           8 = Either space is dense, or MIP is dense
%           10= No progress in this run, return solution from previous one
%           91= Infeasible
%           92= No rectangles to work on
%           93= sLen = 0, no feasible integer solution exists
%
% To make a warm start possible, glcFast saves the following information in
% the file glcFastSave.mat (for internal solver use only):
%   C         Matrix with all rectangle centerpoints, in [0,1]-space.
%   D         Vector with distances from centerpoint to the vertices.
%   F         Vector with function values.
%   G         Matrix with constraint values for each point.
%   Iter      Number of iterations
%   Name      Name of the problem. Used for security if doing warm start
%   Split     Split(i,j) = # splits along dimension i of rectangle j
%   Tr        T(i) is the number of times rectangle i has been trisected.
%   fMinIdx   Indices of the currently best points
%   fMinEQ    sum(abs(infeasibilities)) for minimum points, 0 if no equalities
%   glcfMin   Best function value found at a feasible point.
%   feasible  Flag indicating if a feasible point has been found.
%   ignoreIdx Rectangles to be ignored in the rect. selection procedure.
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
% Direct solver call:
%    Result = glcFast(Prob);
%    PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('glcFast',Prob,2);
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
% To make a restart, just set the restart flag, and call glcFast once again:
%
%    Prob.WarmStart = 1;
%    Note - maximal number of iterations or function evaluations might have
%    to be increased:
%    Prob.optParam.MaxFunc = 20000; % A total of 20000 including previous run
%    Result = tomRun('glcFast',Prob,2);
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
% Recommendation: Put the integers as the first variables !
% Put low range integers before large range integers
% Linear constraints are specially treated
% Equality constraints are added as penalties to the objective.
% Weights are computed automatically, assuimg f(x) scaled to be roughly 1
% at optimum. Otherwise scale f(x)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Feb 15, 1999.   Last modified Feb 22, 2007.

function Result = glcFast(Prob)

if nargin < 1
   error('glcFast needs input structure Prob');
end

solvType=checkType('glc');

Prob=ProbCheck(Prob,'glcFast',solvType);

Prob = iniSolve(Prob,solvType,0,0);

MaxCPU = DefPar(Prob,'MaxCPU',1E300);

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
      error('glcFast: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

isMIP = ~isempty(IntVars);

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

% Safeguard
if isempty(MaxIter) | MaxIter < 0
   MaxIter = max(5000,n*1000);
end
if isempty(MaxFunc) | MaxFunc < 0
   MaxFunc = max(10000,n*2000);
end

if isinf(fGoal) | isempty(fGoal),   fGoal = -1E300; end

Result                 = ResultDef(Prob);
Result.Solver          = 'glcFast';
Result.SolverAlgorithm = 'Constrained DIRECT - Fortran Mex implementation';


%if any(isinf(x_L)) | any(isinf(x_U))
%   disp('glcFast solves box-bounded problems.');
%   disp('Found some bound to be Inf');
%   Result.ExitFlag = 2;
%   Result.ExitText =str2mat('glcFast solves box-bounded problems' ...
%                           ,'Found some bound to be Inf');
%   Result=endSolve(Prob,Result);
%   return
%end

%A   = Prob.A;        % Linear constraint matrix
b_L = Prob.b_L(:);   % Lower bounds, linear constraints
b_U = Prob.b_U(:);   % Upper bounds, linear constraints
%c_L = Prob.c_L(:);   % Lower bounds, nonlinear constraints
%c_U = Prob.c_U(:);   % Upper bounds, nonlinear constraints

% Forgetting factor in alfa * s_0 + (1-alfa)* s_0_new

alfa = 0.99;

% Factor in infeasibility sum for equalities

betaFac = 100;

% Number of function evaluations to use for average, before switching to
% exponential forgetting

AvIter  = 100;


fcn  = Prob.FUNCS.f;
cfcn = Prob.FUNCS.c;
if isempty(cfcn)
   if xnargin(fcn) > 1 
      probPtr = Prob;
   else
      probPtr = -999;
   end
else
   n1 = xnargin(fcn);
   n2 = xnargin(cfcn);
   if (n1 == 2 & n2 == 2) | ( n1 == 3 & n2 == 3) | (n1 == 3 & n2 == 2)
      probPtr = Prob;
   elseif n1 == 1 & n2 == 1
      probPtr = -999;
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
      error('glcFast: Same number of input arguments is needed');
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

AxFeas    = DefPar(GD,'AxFeas',0); 
fcALL     = DefPar(GD,'fcALL',0); 
LinWeight = DefPar(GD,'LinWeight',0.1); 
fEqual    = DefPar(GD,'fEqual',1E-10); 
AxFeas    = min(fcALL,AxFeas);  % AxFeas must be 0 if fcALL is 0

IntPar  = [MaxIter,MaxFunc, AvIter,AxFeas];
RealPar = [fTol, EpsGlob, cTol, xTol, bTol, alfa, LinWeight, fEqual, ...
           betaFac, MaxCPU, fIP];


if Prob.WarmStart
   % Restart with values from previous run.

   Name1 = deblank(Prob.Name);  % Problem name
   load('glcFastSave.mat','Name')
   if strcmp(Name1,Name)
      load('glcFastSave.mat','F','Iter')
      nFunc0 = length(F); 
      Iter0 = Iter; 
      clear F
      if PriLev > 0
         fprintf('\n Restart with %d sampled points from last runs\n',nFunc0);
      end
   else
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files glcFastSave.mat?\n');
      end
      nFunc0 = 0; 
      Iter0 = 0; 
      Name = Name1; 
   end
else
   nFunc0 = 0; 
   Iter0 = 0; 
   Name = deblank(Prob.Name);  % Problem name
end

%IterPrint = 2*(PriLev > 1) + IterPrint;
if PriLev > 1
   IterPrint = PriLev;
else
   IterPrint = max(IterPrint,PriLev);
end

[glcfMin, xMin, Iter, nFunc, convFlag, minSplit, minPnts] = tomsol(18, fcn, ...
   cfcn, Prob.WarmStart, IterPrint, x_L, x_U, full(Prob.A), ...
   b_L, b_U, Prob.c_L, Prob.c_U, ...
   IntVars, IntPar, fGoal, RealPar, probPtr);

try
   save 'glcFastSave.mat' Name -APPEND
catch
   warning('Failed to append variable Name to glcFastSave.mat');
   disp(lasterr);
end

switch convFlag
 case 3
   ExitFlag = 0;
   Inform   = 3;
   if Prob.WarmStart
      Result.ExitText = ['Max iterations. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Maximum number of iterations reached';
   end
case 9
   ExitFlag = 0;
   Inform   = 9;
   if Prob.WarmStart
      Result.ExitText = ['CPU time limit reached. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'CPU time limit reached.';
   end
 
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
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Maximum number of f(x) evaluations reached';
   end
 case 1
   ExitFlag = 0;
   Inform   = 1;
   if Prob.WarmStart
      Result.ExitText = ['Function value less than fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Converged to fGoal';
   end
 case 11
   ExitFlag = 0;
   Inform   = 2;
   if Prob.WarmStart
      Result.ExitText = ['Converged to fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Converged to fGoal';
   end
 case -1
   ExitFlag = 0;
   Inform   = 92;
   if Prob.WarmStart
      Result.ExitText = ['No rectangles to work on. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'No rectangles to work on';
   end
 case -4
   ExitFlag = 8;
   Inform   = 93;
   if Prob.WarmStart
      Result.ExitText = ['No feasible integer solution exists. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'No feasible integer solution exists';
   end
 otherwise
   ExitFlag = 4;
   Inform   = 99;
   Result.ExitText = 'Numerical trouble - EXIT!';
end

if minSplit == -1
   Inform   = 7;
   if Prob.WarmStart
      Result.ExitText = ['Space is dense. ', ...
          'Tried ' num2str(nFunc) ' f(x) values'];
   else
      Result.ExitText = 'Space is dense';
   end
end
if minSplit == -2
   % ExitFlag = 0; Does this only occur if convFlag = -1
   Inform   = 8;
   if Prob.WarmStart
      Result.ExitText = ['Space or MIP is dense. ', ...
          'Tried ' num2str(nFunc) ' f(x) values'];
   else
      Result.ExitText = 'Either space is dense, or MIP is dense';
   end
end
if minSplit < 0
   Result.maxTri   = 0;            % Maximal size of possible triangles
else
   Result.maxTri   = 3.^-minSplit; % Maximal size of possible triangles
end
Result.Iter     = Iter-Iter0;   % Number of iterations this run
Result.FuncEv   = nFunc-nFunc0; % Number of function evaluations this run
if Prob.mNonLin > 0
   Result.ConstrEv = Result.FuncEv;
end

load 'glcFastSave.mat' feasible fMinIdx

%if isempty(minPnts)
if ~feasible
   dBIG = 1E300;
   load 'glcFastSave.mat' F C G
   if Inform < 90 
      Inform = 91; 
      ExitFlag = 7;
      Result.ExitText = ['Tried ' num2str(nFunc) ...
                         ' f(x) values in total, no feasible point'];
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
   %keyboard
   try
      save 'glcFastSave.mat' fMinIdx -APPEND
   catch
      warning('Failed to append variable fMinIdx to glcFastSave.mat');
      disp(lasterr);
   end
      
elseif isempty(minPnts)
   % if minPnts [], no progress in this run, get best point from previous run
   load 'glcFastSave.mat' F C G
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
   Result.c_k = c_k; % Constraint value at x_k;
end
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

% Must set global n_f to see number of function evaluations in the GUI

global n_f
n_f = nFunc-nFunc0;

Result          = endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 010730  jkk  Moved everything to Fortran
% 010815  hkh  Clean up
% 010817  hkh  Avoid calling Tomlab gateway nlp_f and nlp_c, call user
%              routines directly. Return max size of triangles as maxTri
% 011030  hkh  Improve comments. Add Name to mat-file. Check of warm start.
% 011131  hkh  Bug fixes for GUI. Revision.
% 011110  hkh  Fixed errors in comments
% 011112  hkh  Add test on fGoal empty, set as very low then, avoid DLL crash
% 020110  hkh  Change name fMin to glcfMin. Conflict in 5.x with fmin function
% 020311  hkh  Use bTol for linear constraints
% 020312  hkh  Add LinWeight parameter balance linear to nonlinear constraints
% 020325  hkh  Add fEqual parameter, tolerance when points are equal
% 020326  hkh  Skip sending the length of the nonlinear constraints
% 020326  hkh  Add output minPnts with indices to optimal points
% 030124  hkh  Allow general IntVars input similar to mipSolve
% 030124  hkh  Do preSolve if any lower or upper bound is inf or empty
% 030124  hkh  Accept empty lower or upper bound as inpu
% 030129  hkh  Return convFlag in Result.Inform
% 030425  hkh  Change comments and adding more print levels
% 040111  hkh  Change call to inisolve
% 040330  hkh  Revised Inform values, revised ExitFlag values
% 040331  hkh  If infeasible, find point close to feasible
% 040331  hkh  Always return nonlinear constraint values Result.c_k
% 040403  hkh  Return previous solution, if no progress this run
% 040415  hkh  Change ExitFlag from 4 to 8 for Inform=92
% 040428  hkh  Set maxTri=0 if minSplit < 0
% 040517  ango Send MaxCPU to tomsol MEX. Handle return code for time limit
% 041017  hkh  Add input xIP and fIP, Use fIP when computing f_star
% 041018  ango Minor issues fixed related to 041017 change. 
% 041123  hkh  Change call to tomRun in help
% 050211  hkh  Test if # of inputs different in fcn and cfcn
% 050216  hkh  ExitFlag=0 if convFlag=-1 or minSplit=-2
% 050314  ango Add Prob.glcDirect.AxFeas
% 050405  hkh  Move LinWeight, fEqual from Prob.GO to Prob.glcDirect
% 050405  hkh  Add fcALL to Prob.glcDirect, change AxFeas
% 051006  hkh  Safe handling of MaxFunc, new default
% 060327  hkh  Add output of Result.ConstrEv=Result.FuncEv, if mNonLin > 0
% 060327  hkh  Allow both 2 and 3 user function input arguments (3rd=varargin)
% 060814  med  FUNCS used for callbacks instead
% 061003  ango try-catch on save statement
% 061122  ango More try-catch
% 070222  hkh  Revise IntVars handling, use new format
