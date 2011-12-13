% sepMINLP.m
%
% sepMINLP solves general constrained mixed integer nonlinear non-convex
% optimization by a separable approach.
% A nonlinear problem that is a nonlinear least squares type of problem
% is treated with special solvers.
%
% The problem given as input in Prob is separated into one IP problem,
% defined in ProbIP and one MINLP/NLP/NLLS problem, defined in ProbIP.Prob
%
% The first linear constraint should only have 0s and 1s as coefficients.
% The variables in the IP problem is defined by the first linear constraint
% as the variables with coefficients 1, defined as the set IP.
% The lower and upper bound of this constraint defines the number of
% variables in IP that must be assigned to 1.
% The algorithm for problem IP assures that no subproblem is solved when
% the solution cannot be feasible with respect to constraint 1 and other
% linear constraints dependent only on the variables in the set IP.
%
% The trial solutions in the integer part (IP) are generated by either
% the solver glcSolve or the recursive solver genIP.
%
% There are two approaches to solve the MINLP/NLP/NLLS subproblem:
%
% 1.  The subproblem is changed and shrinked dependent on what integers
%     are selected to be 1.
%     The names of the modified functions are given in Prob.MIP.FUNCS
% 2.  The subproblem is the full problem, but with the variables in the set IP
%     set to fixed values (either 0 or 1). Used if Prob.MIP.FUNCS == [].
%
% ----------------
% Solver glcSolve:
% ----------------
%
% The routine sepMIPD_f is called to solve the MINLP/NLP/NLLS problem
% using standard TOMLAB calls, i.e. using nlp_f.m as the TOMLAB gateway.
% Communication from sepMINLP to sepMIPD_f is using the global variables
% allIP, feasIP and fBest.
%
% -------------
% Solver genIP:
% -------------
%
% genIP calls the sub function sepMIP_f to solve the MINLP/NLP/NLLS
% problem.
%
%   sepMIP_f(IntVars)
%
% genIP stops the computation if CPU time used > MaxCPU.
%
% genIP is a recursive generation of all integers with RMmax ones and the
% rest of values 0.
%
% Communication from sepMINLP to sepMIP_f is using the global variables
% allIP, feasIP, fBest, ResBest and ProbIP to send the input in ProbIP
% and pick up the best result obtained from sepMIP_f in fBest and ResBest
% allIP. feasIP is incremented in sepMIP_f.
%
% allIP is the sum of all integer combinations tried.
% feasIP is the sum of all integer combinations that are integer feasible,
% where sepMIP_f has tried to solve the subproblem.
%
% --------------------------------------------------------------------------
% sepMINLP solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Calling syntax:
%
% function Result = sepMINLP(Prob)
%
% INPUT PARAMETERS
%
% Prob    Structure, where the following variables are used:
%
%   Name      Name of the problem. Used for security if doing warm start
%   x_L       Lower bounds for each element in x. Try to set tight bounds.
%             Infinity is replaced by -10000.
%   x_U       Upper bounds for each element in x. Try to set tight bounds.
%             Infinity is replaced by 10000.
%
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix. The first linear constraint is used
%             for the separation
%
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%
%   PriLevOpt Print Level
%
%   MaxCPU    Maximal CPU Time (in seconds) to be used
%
% Solver      Structure with choice of solver
%   Alg       = 1        Use glcSolve
%             otherwise: Use genIP
%
% optParam    Structure in Prob, Prob.optParam.  (USED BY glcSolve)
%             Defines optimization parameters. Fields used:
%   MaxIter   Maximal number of iterations, default max(5000,n*1000);
%   MaxFunc   Maximal number of function evaluations, default max(10000,n*2000)
%   IterPrint Print one line each iteration
%   cTol      Nonlinear constraint feasibility tolerance
%   bTol      Linear constraint feasibility tolerance
%   fGoal     Goal for function value, if empty not used
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal ~=0
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
%             If any element >1, IntVars is the indices for integer
%             variables.
%
%   fIP       An upper bound on the optimal f(x) value. If empty, set as Inf.
%   xIP       The x-values giving the fIP value.
%             If fIP empty and xIP given, fIP will be computed
%             if xIP nonempty, its feasibility is checked
%
%   FUNCS     Subfields similar to Prob.FUNCS. Give the function names
%             for the functions written to handle the SHRINKED subproblem
%             with.
%
% ---------------------------------------
% glcDirect   Structure with DIRECT algorithm specific parameters. Fields used:
% ---------------------------------------
% See help glcDirect for a description of all special parameters used.
%
% OUTPUT PARAMETERS
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
% USAGE:
%
%      Prob.Solver.Alg = 1;
%      Prob.MaxCPU = 100;
%      Result = sepMINLP(Prob);

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 2004-2007 by Tomlab Optimization Inc., $Release: 5.9.0$
% Written Apr 10, 2004.   Last modified May 17, 2007.

function Result = sepMINLP(Prob)

global ResBest IntBest ProbIP;

if nargin < 1
   error('sepMINLP needs input structure Prob');
end

solvType=checkType('glc');

Prob=ProbCheck(Prob,'glcSolve',solvType);

Prob = iniSolve(Prob,solvType,0,0);

% Pick up input parameters from the Prob structure:
x_L = Prob.x_L(:);   % Lower bounds
x_U = Prob.x_U(:);   % Upper bounds

if isempty(x_L) | isempty(x_U) | any(isinf(x_L) | isinf(x_U))
   Prob2 = preSolve(Prob);
   x_L = Prob2.x_L;
   x_U = Prob2.x_U;
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
if isempty(IntVars)
   error('Empty set of integer variables. sepMINLP only solves MINLP');
end
% Always IntVars as index vector, for use in subproblem
Prob.MIP.IntVars = IntVars;
% Field for shrinked functions
if ~isfield(Prob.MIP,'FUNCS')
   % Define empty field for use in sepMIP_f
   Prob.MIP.FUNCS  = [];
end

% FUNCS = Prob.MIP.FUNCS;
% PriLev    = Prob.PriLevOpt;          % Print level
% IterPrint = Prob.optParam.IterPrint; % Print short information each iteration

% fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
% if isinf(fGoal),   fGoal = -1E300; end
% if isempty(fGoal), fGoal = -1E300; end

if isfield(Prob.Solver,'Alg')
   Solver   = Prob.Solver.Alg;
else
   Solver = [];
end
if isempty(Solver), Solver = 0; end

Result = ResultDef(Prob);

% Basic ProbIP holds the IP problem, ProbIP.Prob holds the NLP problem
ProbIP       = Prob;
% Determine the integer variable set IP for problem IP
IP           = find(Prob.A(1,:) == 1);
nI           = length(IP);
if nI == 0
   error('Illegal 1st linear constraint. Empty set of integer variables');
end
ProbIP.N      = nI;
ProbIP.x_L    = x_L(IP);
ProbIP.x_U    = x_U(IP);
ProbIP.x_0    = x_L(IP);
ProbIP.c_L    = [];
ProbIP.c_U    = [];
ProbIP.x_max  = [];
ProbIP.x_opt  = [];
ProbIP.f_opt  = [];
ProbIP.MIP.IP = IP;
% Logical vector for variables not in IP
IV           = true(n,1);
IV(IP)       = 0;
% ix - index vector for variables not in IP
ix           = find(IV);
% i0 - logical vector, 1 = linear constraints independent of non-IP vars
i0           = full(all(Prob.A(:,ix)' == 0));

if i0(1) == 0
   error(...
    'Illegal 1st linear constraint. Nonzero element for non-integer variable');
end
if length(i0) > 1
   % iv - index for additional linear constraints only dependent of IP vars
   iv               = 1+find(i0(2:end));
   ProbIP.MIP.A     = Prob.A(iv,IP);
   ProbIP.MIP.b_L   = Prob.b_L(iv);
   ProbIP.MIP.b_U   = Prob.b_U(iv);
else
   iv               = [];
   ProbIP.MIP.A     = [];
   ProbIP.MIP.b_L   = [];
   ProbIP.MIP.b_U   = [];
end
%A  = full(ProbIP.MIP.A)
j  = length(iv);
% Check for constraints with 1 element in non-separable part
% i1 - index for linear constraints with 1 element in non-IP, the rest in IP
i1 = find(full(sum(Prob.A(:,ix)' ~= 0) == 1 & any(Prob.A(:,IP)' ~= 0)));
nA = length(i1);
if nA > 0
   % Add constraints with corrected lower/upper bound
   % If x has tight bounds, these can avoid some unnecessary calls
   [ii,jj] = find(Prob.A(i1,ix)');
   jx   = ix(ii);
   % [i1(:),jx]
   ProbIP.MIP.A(j+1:j+nA,:) = Prob.A(i1,IP);
   for i = 1:nA
       i1i = i1(i);
       jxi = jx(i);
       ProbIP.MIP.b_L(j+i) = ...
       min(Prob.b_L(i1i)-Prob.A(i1i,jxi)*x_L(jxi),...
           Prob.b_L(i1i)-Prob.A(i1i,jxi)*x_U(jxi));
       ProbIP.MIP.b_U(j+i) = ...
       max(Prob.b_U(i1i)-Prob.A(i1i,jxi)*x_L(jxi),...
           Prob.b_U(i1i)-Prob.A(i1i,jxi)*x_U(jxi));
   end
   % xprint([ProbIP.MIP.b_L(j+1:end) jx ProbIP.MIP.b_U(j+1:end)]',[],[],3 )
else
   jx = [];
end
Prob.CHECK=0;
if length(IntVars) == nI
   probType          = checkType('con');
   ProbIP.Solver     = 'snopt';
   Prob.optParam     = optParamDef('snopt',probType);
   Prob              = ProbCheck(Prob,'snopt',probType);
   ProbIP.Prob.Name  = [Prob.Name ' - NLP part'];
else
   probType = checkType('minlp');
   ProbIP.Solver     = 'minlpBB';
   Prob.optParam     = optParamDef('minlpBB',probType);
   Prob              = ProbCheck(Prob,'minlpBB',probType);
   ProbIP.Prob.Name  = [Prob.Name ' - MINLP subproblem'];
end
ProbIP.Prob       = Prob;
ProbIP.Name       = [Prob.Name ' - separable IP part'];
nLow              = Prob.b_L(1);
nUpp              = Prob.b_U(1);
if isempty(Prob.MIP.FUNCS)
else
   % Linear constraints involving nonseparated variables
   % with more than one nonzero coefficient
   i2L            = ones(length(i0),1);
   i2L(1)         = 0;
   i2L(iv)        = 0;
   i2L(i1)        = 0;
   i2             = find(i2L);
   if ~isempty(i2)
      % Additional linear constraints added in special field
      ProbIP.Prob.LINCON.A   = Prob.A(i2,ix);
      ProbIP.Prob.LINCON.b_L = Prob.b_L(i2);
      ProbIP.Prob.LINCON.b_U = Prob.b_U(i2);
   else
      ProbIP.Prob.LINCON  = [];
   end
   ProbIP.MIP.i1 = [i1(:),jx];
   ProbIP.MIP.i2 = i2(:);
end

ProbIP.MIP.IntVars = 1:nI;
ProbIP.mNonLin     = 0;

global allIP feasIP fBest
feasIP  = 0;
allIP   = 0;
if isempty(Prob.MIP.fIP)
   fBest = inf;
else
   fBest = Prob.MIP.fIP;
end

Prob.SaveRes = 1;
if Solver == 1
   ProbIP = tomFiles(ProbIP, 'sepMIPD_f', [], [], [], [], [], [], [], [], [], []);
   AA = ProbIP.MIP.A;
   if isempty(AA)
      ProbIP.A=ProbIP.A(1,:);
      ProbIP.b_L=ProbIP.b_L(1);
      ProbIP.b_U=ProbIP.b_U(1);
   else
      % Add logical constraints
      ProbIP.A=[ProbIP.A(1,:);AA];
      ProbIP.b_L=[ProbIP.b_L(1);ProbIP.MIP.b_L];
      ProbIP.b_U=[ProbIP.b_U(1);ProbIP.MIP.b_U];
   end
   ProbIP.mLin = size(ProbIP.A,1);
   Result=tomRun('glcSolve',ProbIP,[],2);
   %rsnopt( [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
   %            [], [], [], [], [], [], [], 5 );
   % Result         = snoptTL(Result.SubResult.Prob,5);
   Result.Solver          = 'sepMINLP';
   %Result.f_k             = Result.SubResult.f_k;
   %Result.c_k             = Result.SubResult.c_k;
   Result.SolverAlgorithm = 'Separable MINLP: glcSolve for IP+Local NLP Solver';
elseif Solver == 2
   % Compute logical constraints in mipd_f2.m
   ProbIP = tomFiles(ProbIP, 'mipd_f2', [], [], [], [], [], [], [], [], [], []);
   ProbIP.A=ProbIP.A(1,:);
   ProbIP.b_L=ProbIP.b_L(1);
   ProbIP.b_U=ProbIP.b_U(1);
   ProbIP.mLin = size(ProbIP.A,1);
   Result=tomRun('glcSolve',ProbIP,[],2);
   Result.Solver          = 'sepMINLP';
   Result.SolverAlgorithm = 'Separable MINLP: glcSolve IPII+Local NLP Solver';
   Result.IPIter          = feasIP;      % Number of iterations this run
   Result.IPTotal         = allIP;       % Total iterations, incl.skipped
else
   % Solve IP with genIP, subfunction sepMIP_f
   MaxCPU                 = Prob.MaxCPU;
   TIME00                 = Prob.TIME0;
   ResBest                = Result;
   IntBest                = [];
   endTime                = TIME00+MaxCPU;
   for i = nUpp:-1:nLow
       genIP('',nI,i,endTime)
   end
   Result                 = ResBest;
   Result.IPIter          = feasIP;      % Number of iterations this run
   Result.IPTotal         = allIP;       % Total iterations, incl.skipped
   Result.IntVars         = IntBest;
   Result.Solver          = 'sepMINLP';
   if length(IntVars) == nI
      Result.SolverAlgorithm = 'Separable MINLP: genIP for IP+Local NLP Solver';
   else
      Result.SolverAlgorithm = 'Separable MINLP: genIP for IP+MINLP subproblem';
   end
   if ~isinf(fBest) 
      %ProbIP.Prob          = Result.Prob;
   else
      ProbIP.Prob          = [];
   end
end

Result = endSolve(ProbIP,Result);

function genIP(s,m,n,endTime)

if cputime > endTime 
   fprintf('genIP: Time limit reached\n');
   return; 
end
if m==0 
   sepMIP_f(find(s=='1'))
elseif n==0
   sepMIP_f(find(s=='1'))
else
   genIP( [ s '1' ] , m-1, n-1, endTime );
   if m>n ,genIP( [ s '0' ] , m-1, n, endTime   ); end
end

% MODIFICATION LOG:
%
% 040426  hkh  Written
% 070517  hkh  Generalized to general separable MINLP solver