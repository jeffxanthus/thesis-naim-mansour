% FilMINT.m:
%
% function Result = FilMINT(Prob)
%
% An Outer Approximation-Based solver for 
% Mixed-Integer Nonlinear Programming (MINLP). 
% Handle both convex or nonconvex sub problems.
%
% The parameter Convex (see below) determines if to assume the NLP subproblems are
% convex or not.
%
% Solving the following MINLP 
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds
%    b_L <= A*x  <= b_U, linear constraints
%    c_L <= c(x) <= c_L, nonlinear constraints
%
%    x_i are restricted to integer values for i in I,
%        where I is a subset of {1,2,...,n}, x in R^n
%
% INPUT PARAMETERS
% Fields in Prob:
%   x_L:    Lower bounds on x. If empty assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%   A:      The linear constraint matrix 
%   b_L:    Lower bounds on linear constraints. 
%   b_U:    Upper bounds on linear constraints
%   c_L     Lower bounds on nonlinear constraints
%   c_U     Upper bounds on nonlinear constraints
%
%   x_0:    Starting point x 
%
% Convex    If Convex==1, assume NLP problems are convex, and only one local
%           NLP solver call is used at each node.
%           If Convex==0, multiMin is used to do many calls to a local solver to
%           determine the global minima at each node. The global minimum with
%           most components integer valued is chosen. (DEFAULT)
%
% MaxCPU:   Maximal CPU Time (in seconds) to be used by FilMINT,
%           stops with best point found
% PriLev    Print level in FilMINT (default 1). Also see optParam.IterPrint
% PriLevOpt Print level in sub solvers (SNOPT and other NLP solvers):
%           =0 No output; >0 Convergence results; 
%           >1 Output every iteration  >2 Output each step in the NLP alg
%           For other NLP solvers, see the documentation for the solver
% HKH Should be print in LP solver as well? Maybe not.
%
% WarmStart If true, >0, FilMINT reads the output from the last
%           run from Prob.FilMINT, if it exists. If it doesn't exist, FilMINT attempts 
%           to open and read warm start data from mat-file FilMINTSave.mat. 
%           FilMINT uses the warm start information to continue from the last run.
%           The mat-file FilMINTSave.mat is saved every Prob.MIP.SaveFreq iteration
%
% SolverNLP Name of the solver used for NLP subproblems. If empty,
%           the default solver is found calling GetSolver('con',1);
%           If TOMLAB /SOL installed, SNOPT is the default solver.
%           If SolverNLP is a SOL solver (SNOPT, MINOS or NPSOL), 
%           the SOL.optPar and SOL.PrintFile is used:
%           See help minosTL.m, npsolTL.m or snoptTL.m for how to set these parameters
%
% SolverLP  Name of the solver used for LP subproblems. If empty,
%           the default solver is found calling GetSolver('lp',1);
%           If TOMLAB /SOL or /SNOPT or /NPSOL installed, MINOS is the default solver.
%           If SolverLP is a SOL solver (MINOS, SNOPT or NPSOL), 
%           the SOL.optPar and SOL.PrintFile is used:
%           See help minosTL.m, npsolTL.m or snoptTL.m for how to set these parameters
% HKH NOTE - If to change LP settings with SOL.optPar, may be conflict with NLP settings
%           Maybe use other field instead, if necessary to change LP settings, maybe Prob.MIP.SOL.optPar?
% HKH - see mipSolve for how to use MINOS efficiently
%
% RandState If Convex == 0, RandState is sent to multiMin to initialize the random generator
%           RandState is used as follows:
%           If > 0, rand('state',RandState) is set to initialize the pseudo-random generator
%           if < 0, rand('state',sum(100*clock)) is set to give a new set
%           of random values each run
%           if RandState == 0, rand('state',) is not called
%           Default RandState = -1
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%           Defines integer optimization parameters. Fields used:
% IntVars:  
%           If empty, all variables are assumed non-integer 
%           If islogical(IntVars) (=all elements are 0/1), then
%           1 = integer variable, 0 = continuous variable.
%           If any element >1, IntVars is the indices for integer variables
%
% VarWeight:Weight for each variable in the variable selection phase.
%           A lower value gives higher priority. Setting
%           Prob.MIP.VarWeight might improve convergence.
% DualGap   FilMINT stops if the duality gap is less than DualGap
%           DualGap = 1, stop at first integer solution
%           e.g. DualGap = 0.01, stop if solution < 1% from optimal solution
% fIP       An upper bound on the IP value wanted. Makes it possible to
%           cut branches and avoid node computations. Used even if xIP not given
% xIP       The x-values giving the fIP value, if a solution (xIP,fIP) is known.
%
% NodeSel   Node selection method in branch and bound
%           = 0 Depth First. Priority on  nodes with more integer components.
%           = 1 Breadth First. Priority on  nodes with more integer components.
%           = 2 Depth First. When integer solution found, use NodeSel = 1 (default)
%           = 3 Pure LIFO (Last in, first out) Depth First
%           = 4 Pure FIFO (First in, first out) Breadth First
%           = 5 Pure LIFO Depth First. When integer solution found, use NodeSel 4
% HKH implement more NodeSel strategies:
%   -    Best bound
%   -    Best estimate
%   -    Adaptive
% VarSel    Variable selection method in branch and bound
%           = 1 Use variable with most fractional value
%           = 2 Use gradient and distance to nearest integer value 
% HKH Implement more variable selection strategies, 1st: Pseudo cost, later strong branching
% KNAPSACK  If = 1, use a knapsack heuristic. Default 0.
% ROUNDH    If = 1, use a rounding heuristic. Default 0.
%
% SaveFreq  Warm start info saved on FilMINTSave.mat every SaveFreq iteration
%           (default -1, i.e. no warm start info is saved)
% -------------------------------------------------------------------------------
% optParam  Structure in Prob. Fields used in Prob.optParam, also in sub solvers:
% -------------------------------------------------------------------------------
%  MaxIter   Maximal number of iterations, default 10000
%  IterPrint Print short information each iteration (PriLev > 0 ==> IterPrint = 1)
%            Iteration number: Depth in tree (symbol L[] - empty list, 
%            symbol Gap - Dual Gap convergence)
%            fNLP (Optimal f(x) current node), fIPMin (Best integer feasible f(x) found)
%            LowBnd (Lower bound on optimal integer feasible f(x)),
%            Dual Gap in absolut value and percent.
%            The length of the node list L, |L|
%            The Inform and ExitFlag the solver returned at the current node
%            FuEv (Number of function evaluations used by solver at current node)
%            date/time stamp 
%  bTol      Linear constraint violation convergence tolerance
%  cTol      Constraint violation convergence tolerance
%
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag 
%            == 0  => Global optimal solution found, or integer solution with
%                     duality gap less than user tolerance
%            == 1  => Maximal number of iterations reached.
%            == 2  => Empty feasible set, no integer solution found.
%            == 3  => Rank problems (not used)
%            == 4  => No feasible point found running NLP relaxation
%            == 5  => Illegal x_0 found in NLP relaxation ???
%            == 99 => Maximal CPU Time used (cputime > Prob.MaxCPU)
%   Inform   Code telling type of convergence, returned from subsolver
%   ExitText Text string giving ExitFlag and Inform information
%   DualGap  Relative duality gap, max(0,fIPMin-fLB)/|fIPMin|, if fIPMin ~=0;
%            max(0,fIPMin-fLB) if fIPMin == 0. If fIPMin ~=0:
%            Scale with 100, 100*DualGap, to get the percentage duality gap.
%            For absolute value duality gap: scale with fIPMin, fIPMin * DualGap
%   x_k      Solution
%   v_k      Lagrange multipliers. Bounds, Linear and Nonlinear Constraints, n+mLin+mNonLin
%   f_k      Function value at optimum
%   g_k      Gradient vector at optimum
%   x_0      Starting point x_0
%   f_0      Function value at start
%   c_k      Constraint values at optimum
%   cJac     Constraint derivative values at optimum
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%   bState   Linear constraint: Inactive==0; On lower bound == 1; 
%            On upper bound == 2; Equality == 3;
%   cState   State of each general constraint.
%   Solver   FilMINT
%   SolverAlgorithm  Description of method used
%
%   MINLP 
%       A structure with warm start information. Use with WarmDefGLOBAL, see example below
%
% -- Calling FilMINT
%
% Recommended is to first set IterPrint, to get information each iteration
%      Prob.optParam.IterPrint = 1;
%
% Driver call, including printing with level 2:
%      Result = tomRun('FilMINT',Prob,2);
%
% Direct solver call:
%      Result = FilMINT(Prob);
%      PrintResult(Result);
%
% -- Warm start
% To make a restart (warm start), just set the warm start flag, and call
% FilMINT once again:
%
%      Prob.WarmStart = 1;
%      Result         = tomRun('FilMINT', Prob, 2);
%
% FilMINT will read warm start information from the FilMINTSave.mat file
%
% Another warm start (with same MaxFunc) is made by just calling tomRun again
%      Result = tomRun('FilMINT', Prob, 2);
%
% To make a restart from the warm start information in the Result
% structure, make a call to WarmDefGLOBAL before calling FilMINT.
% WarmDefGLOBAL moves information from the Result structure
% to the Prob structure and sets the warm start flag, Prob.WarmStart = 1; 
%
%      Prob = WarmDefGLOBAL('FilMINT', Prob, Result);
% where Result is the result structure returned by the previous run.
% A warm start (with same MaxIter) is done by just calling tomRun again
%      Result = tomRun('FilMINT', Prob, 2);
%
% To make another warm start with new MaxIter 100, say, redefine MaxIter as:
%      Prob.optParam.MaxIter = 100;
% Then repeat the two lines:
%      Prob = WarmDefGLOBAL('FilMINT', Prob, Result);
%      Result = tomRun('FilMINT', Prob, 2);
%


% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written Nov 18, 2009.   Last modified Nov 25, 2009.

function ResultMIP = FilMINT(Prob)

if nargin < 1
   error('FilMINT needs input structure Prob');
end

%HKH2
global NLP_x NLP_f NARG NLP_xc NLP_c NLP_xg NLP_g NLP_xdc NLP_dc 
% Insert reset of global variables before single call to nlp-routines
% if tomRun call with other Prob struct has been done
% Before nlp_f
      % NARG = []; NLP_x=[]; NLP_f=[]; 
% Before nlp_g:
      % NARG = []; NLP_x=[]; NLP_f=[]; NLP_xg=[]; NLP_g=[]; 
% Before nlp_c:
      % NARG = []; NLP_xc=[]; NLP_c=[]; 
% Before nlp_dc:
      % NARG = []; NLP_xc=[]; NLP_c=[]; NLP_xdc=[]; NLP_dc=[]; 

solvType  = checkType('minlp');
Prob      = ProbCheck(Prob,'FilMINT',solvType);
Prob      = iniSolve(Prob,solvType,1,1); % DEPENDS ON SUBSOLVER, SOL ASSUMED
absviol   = 0;   % Check relative error in linear and nonlinear constraints

% Avoid call to estimate Lagrange multipliers in PrintResult
Prob.PrintLM = 0;

Name      = deblank(Prob.Name);
n         = Prob.N;
x_L       = Prob.x_L;
x_U       = Prob.x_U;
x_L0      = x_L;
x_U0      = x_U;
b_L       = Prob.b_L;
b_U       = Prob.b_U;
c_L       = Prob.c_L;
c_U       = Prob.c_U;
A         = Prob.A;
dLin      = Prob.mLin;
dNoL      = Prob.mNonLin;
m         = dLin + dNoL;
nb        = n + dLin + dNoL;
%HKH Should check if MILP and then set Convex = 1; directly
% if FUNCS.f is lp_f or qp_f or lls_f and Prob.mNonLin == 0
% User guess if the problem is convex or not
CONVEX    = DefPar(Prob,'Convex',0);
% Limit on CPU time
MaxCPU    = DefPar(Prob,'MaxCPU',Inf);

optParam  = Prob.optParam;
MaxIter   = optParam.MaxIter;  % Maximal number of iterations
wait      = optParam.wait;     % Pause after printout if true, when DEBUG==1
bTol      = optParam.bTol;     % Linear constraint violation convergence tolerance
cTol      = optParam.cTol;     % Constraint violation convergence tolerance
% CHeck later on worst case for MxAlloc, maybe to be lowered a bit
MxAlloc   = MaxIter*2+2;

PriLev    = Prob.PriLev;          % Print level
PriLevOpt = Prob.PriLevOpt;       % Print level in subsolvers
IterPrint = optParam.IterPrint;   % Print short information each iteration

if isempty(PriLev), PriLev    = 0; end
if PriLev > 0,      IterPrint = 1; end

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);
% Logical vector for integers
IV       = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('FilMINT: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
% RealVars = find(~IV);
IntVars  = find(IV);
nI       = length(IntVars);
eps_I    = 1E-8;  % Level when a variable is considered as an integer value

ResultMIP                 = ResultDef(Prob);
ResultMIP.Solver          = 'FilMINT';
ResultMIP.SolverAlgorithm = 'Outer Approximation MINLP';
ResultMIP.MINLP           = [];

% fIPMin is upper bound on the optimal solution
fIP    = DefPar(Prob.MIP,'fIP',[]);
fIPMin = fIP;
xIPMin = DefPar(Prob.MIP,'xIP',[]);
if isempty(fIPMin)
   fIPMin = Inf;     % Current BEST integer solution
   if ~isempty(xIPMin)
      % User has supplied an x
      if length(xIPMin) < n
         xIPMin=[xIPMin;zeros(n-length(xIPMin))];
      end
      xIPMin(IntVars) = round(xIPMin(IntVars));
      % Compute the function value at the given point
      fIPMin = nlp_f(xIPMin, Prob);
   end
else
   if isempty(xIPMin)
      fIPMin = Inf;     % Current BEST integer solution
   else
      xIPMin(IntVars) = round(xIPMin(IntVars));
   end
end
% Check if xIPMin is feasible
if ~isinf(fIPMin) 
   if dLin > 0
      ResultMIP.Ax = A*xIPMin;
   end
   if dNoL > 0
       cx          = nlp_c(xIPMin,Prob);
   else
       cx          = [];
   end
   ResultMIP.c_k = cx;
   ResultMIP.x_k = xIPMin;
   % HKH cErrCompute returns || ||_inf
   ResultMIP     = cErrCompute(ResultMIP);
   if any(xIPMin < x_L | xIPMin > x_U) | ...
      ResultMIP.bErr > bTol | ResultMIP.cErr > cTol 
      if IterPrint 
         fprintf('Max error in linear constraints    %f\n',ResultMIP.bErr);
         fprintf('Max error in nonlinear constraints %f\n',ResultMIP.cErr);
         fprintf('xIPMin is not feasible\n');
      end
      fIPMin = Inf;     
      xIPMin = [];     
   end
end
if ~isempty(fIPMin) & isempty(xIPMin)
   fIPCut = fIP;
else
   fIPCut = Inf;
end

if isfield(Prob.MIP,'DualGap')
   if isempty(Prob.MIP.DualGap)
      DualGap = 0;
   else
      DualGap = Prob.MIP.DualGap;
   end
else
   DualGap = 0;
end

VarWeight = DefPar(Prob.MIP,'VarWeight',[]);


% Warm start info saving frequency
SaveFreq           = DefPar(Prob.MIP,'SaveFreq',-1);

% Node selection method
NodeSel            = DefPar(Prob.MIP,'NodeSel',2);

% Variable selection method
VarSel            = DefPar(Prob.MIP,'VarSel',1);

% Heuristics
KNAPSACK          = DefPar(Prob.MIP,'KNAPSACK',0);
ROUNDH            = DefPar(Prob.MIP,'ROUNDH',0);

Prob.MIP.NodeSel  = NodeSel;
Prob.MIP.VarSel   = VarSel;
Prob.MIP.KNAPSACK = KNAPSACK;
Prob.MIP.ROUNDH   = ROUNDH;

if KNAPSACK > 0
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' Knapsack heuristic.'];
   % Check if problem is on standard form or not
   bMK = find(b_L==b_U);
   nMK = length(bMK);
   if nMK > 0  % Find slack variables in equality constraints
      xMK        = zeros(nMK,1);
      [ix jx]    = find(sparse(Prob.A));
      for j = 1:nMK
          k      = bMK(j);
          xMK(j) = max(jx(ix==k)); % Assume slacks are last in row
      end
      xIK = intersect(IntVars,xMK);
   end
end


switch NodeSel
case 0
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Depth First.'];
case 1
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Breadth First.'];
case 3
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' LIFO Depth First.'];
case 4
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' FIFO Breadth First.'];
case 5
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' LIFO Depth First, then FIFO Breadth.'];
otherwise
   NodeSel = 2;
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' Depth First, then Breadth.'];
end

% Variable priority weights
if ~isempty(VarWeight)
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' Priority weights.'];
end

ResultMIP.f_0  = [];
ResultMIP.x_k  = [];
ResultMIP.f_k  = [];
ResultMIP.v_k  = [];

Prob.SolverNLP = DefPar(Prob,'SolverNLP',[]);
if isempty(Prob.SolverNLP)
   SolverNLP=GetSolver('con',1);
   Prob.SolverNLP = SolverNLP;
else
   SolverNLP=Prob.SolverNLP;
end

% HKH Check that MINOS comes first if any SOL toolbox
Prob.SolverLP = DefPar(Prob,'SolverLP',[]);
Prob.SolverLP = 'minos';
if isempty(Prob.SolverLP)
   SolverLP=GetSolver('lp',1);
   Prob.SolverLP = SolverLP;
else
   SolverLP=Prob.SolverLP;
end

% Setup Structure used in NLP calls
ProbNLP             = Prob; % SAME FOR NOW, IntVars ignored.
ProbNLP.MIP.IntVars = [];
ProbNLP.b_L         = b_L; 
ProbNLP.b_U         = b_U;
ProbNLP.A           = A;
ProbNLP.x_L         = x_L; 
ProbNLP.x_U         = x_U; 
ProbNLP.WarmStart   = 0;    % Important to avoid warm start in subsolver

% HKH Setup LP Structure used in NLP calls
ProbLP.c            = [zeros(1,n) 1];
ProbLP.A            = zeros(m+1,n+1);
%HKH2 Must be bug here, m should be n - gives crash on problem 3
%ProbLP.A(1,m+1)     = -1;
ProbLP.A(1,n+1)     = -1;
ProbLP.b_L          = -Inf*ones(m+1,1);
ProbLP.b_U          = zeros(m+1,1);
ProbLP.x_L          = [x_L; -Inf];
ProbLP.x_U          = [x_U;  Inf];
ProbLP.N            = n+1;
ProbLP = lpAssign(ProbLP.c, ProbLP.A, ProbLP.b_L, ProbLP.b_U, ...
                  ProbLP.x_L, ProbLP.x_U, [], 'FilMINT CMP');
% Values for f,c,g and J inserted according to the pattern shown below 
% for each node before solving the continuous master problem (CMP).
% c  = [ 0..0 1]      1  x n+1 
% A  = [ g'  -1      m+1 x n+1
%        J    0 ]P
% bU = [ g'x_k-f_k ; m+1 x 1 
%        J*x_k-c_k ]
% x_L and x_U also need to be updated.

% FilMINT warm start parameters
% HKH Get the warm start structure here - name easily changed
MINLP = DefPar(Prob,'MINLP',[]);

if Prob.WarmStart == 1                    % Restart with values from previous run.

   Name1 = Name;  % Problem name
   
   % If the WarmStart structure exists, then use that warm start data.
   if ~isempty(MINLP)
      Name = MINLP.Name;
      if strcmp(Name1,Name)
         if PriLev >= 2
            fprintf(['Restart with %d iterations from last ' ...
                  'runs.\nUsing warm start data from structure: ' ...
                  'Prob.MINLP.\n'],...
                  MINLP.TotIter);
         end
      else
         Prob.WarmStart = 0;
         if PriLev >= -1000
            fprintf('Previous run was with Problem %s\n',Name);
            fprintf('This run is with Problem %s\n',Name1);
            fprintf(['Impossible to do restart using data from ' ...
                     'structure: Prob.MINLP.\n']);
         end
      end

   else % If no WarmStart structure exists, look for a file.
      filename = 'FilMINTSave.mat';
      if exist(filename,'file')~=2
         fprintf(['Couldn''t find warm start info structure nor warm ' ...
                  'start file in path.\nNo warm start possible.\n']);
         fprintf('exist(%s) gives %d\n',filename,exist(filename))
         Prob.WarmStart = 0;
      else
         load(filename,'Name');

         if strcmp(Name1,Name)
            MINLP = load(filename,'xIPMin','fIPMin','fIPIter','vIPMin',...
                         'FuncEv','GradEv','HessEv', 'ConstrEv', ...
                         'XX','HS','nS', ...
                         'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                         'Depth','Icomp',...
                         'Gap','SumMinIt','RandState','time');
        
            if PriLev >= 2
               fprintf(['Restart with %d iterations from last runs.\n',...
                        'Using warm start data from file: %s.\n'],MINLP.TotIter, filename);
            end
         
         else
            Prob.WarmStart = 0;
            if PriLev >= -1000
               fprintf('Previous run was with Problem %s\n',Name);
               fprintf('This run is with Problem %s\n',Name1);
               fprintf('Impossible to do restart using data from file: %s.\n', filename);
            end
        end
     end
   end
end

ResultMIP.Prob      = Prob;
WarmStart           = Prob.WarmStart;

if WarmStart == 1 
   % Get information from last run
   xIPMin    = MINLP.xIPMin;
   fIPMin    = MINLP.fIPMin;
   fIPIter   = MINLP.fIPIter;
   vIPMin    = MINLP.vIPMin;
   FuncEv    = MINLP.FuncEv;
   GradEv    = MINLP.GradEv;
   HessEv    = MINLP.HessEv;
   ConstrEv  = MINLP.ConstrEv;
   L         = MINLP.L;
   NODE      = MINLP.NODE;
   TotIter   = MINLP.TotIter;
   % Gap       = MINLP.Gap; % Not needed, recomputed below
   SumMinIt  = MINLP.SumMinIt;
   RandState = MINLP.RandState;
   time      = MINLP.time;
   ProbNLP.RandState = RandState;
   ProbNLP.xInit     = max(100,10*ProbNLP.N);
   if isempty(MINLP.HS)
      SOL=0;
   elseif CONVEX
      if strcmpi(SolverNLP,'minos') | strcmpi(SolverNLP,'sqopt') | ... 
          strcmpi(SolverNLP,'lssol') | strcmpi(SolverNLP,'nlssol') | ...
          strcmpi(SolverNLP,'npsol') | strcmpi(SolverNLP,'snopt') | ...
          strcmpi(SolverNLP,'qpopt') 
         SOL=1;
         ProbNLP.WarmStart=1;
      end
   else
      SOL=0;
   end
   % HKH - IMPORTANT
   % Should we pack the tree at warm start or not.
   % In mipSolve switched back to non-packing, but what to do here?????

   % k is the number of active nodes in the list L
   k                  = length(L);
   % Generate the nodes still needed for the active list of nodes
   pred               = MINLP.pred(1:NODE,1);
   p                  = 1;   % Save root node as first node in the tree                     
   % Generate the search tree for each active node
   for i=1:k
       j = L(i);
       while pred(j) > 1     % Generate nodes preceding list node L(i)
          j = pred(j);
          p = [j;p];
       end
   end
   % Find the unique nodes to keep, sorted, will be nodes 1,...,M in new numbering
   u             = unique([p',L]);
   M             = length(u);
   % Allocate arrays
   if SOL > 0
      XX              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      HS              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      nS              = zeros(1,M+MxAlloc);
      XX(:,1:M)       = MINLP.XX(:,u);
      HS(:,1:M)       = MINLP.HS(:,u);
      nS(1,1:M)       = MINLP.nS(:,u);
   else
      XX              = spalloc(n,M+MxAlloc,min(40000,n*MaxIter));
      XX(:,1:M)       = MINLP.XX(:,u);
      HS              = [];
      nS              = [];
   end

   % Lowest value possible at node, from NLP relaxation
   f_min                = Inf*ones(M+MxAlloc,1);      
   f_min(1:M,1)         = MINLP.f_min(u,1);
   % Depth in tree for each node
   Depth                = zeros(M+MxAlloc,1);      
   Depth(1:M,1)         = MINLP.Depth(u,1);
   % Number of integer components
   Icomp                = zeros(M+MxAlloc,1);      
   Icomp(1:M,1)         = MINLP.Icomp(u,1);

   pred                 = zeros(M+MxAlloc,1);
   Problem.xVar         = zeros(M+MxAlloc,1);
   Problem.Bound        = zeros(M+MxAlloc,1);
   Problem.xBound       = zeros(M+MxAlloc,1);

   Problem.xVar(1:M,1)  = MINLP.Problem.xVar(u,1);
   Problem.Bound(1:M,1) = MINLP.Problem.Bound(u,1);
   Problem.xBound(1:M,1)= MINLP.Problem.xBound(u,1);
   % New node numbers Lnew for active nodes in list L
   [Lval,iu,Lnew]       = intersect(L,u);
   % Find new node number for each preceding node
   p                    = MINLP.pred(u,1);
   pred(1)              = 0;                    % Root node
   for i=2:length(p)
       pred(i)          = find(p(i)==u);        % New node numbers in pred
   end
   L                    = Lnew;                 % New node numbers for active list
   NODE                 = M;                    % The number of nodes after reduction
   % Update fLB, lower bound on the solution
   fLB                  = min(f_min(L));

   Inform               = 0;
   ExitFlag             = 0;

   Iter                 = 0;
   SumIter              = TotIter;              % SumIter + Iter always total iterations

else                                            % end of warm start
   % Initial step
   WarmStart            = 0;
   if SaveFreq >= 0
      if exist('FilMINTSave.mat','file')
         delete('FilMINTSave.mat');
      end
   end

   if CONVEX
      ResultNLP = tomRun(SolverNLP,ProbNLP,PriLev > 0);
      ExitFlag = ResultNLP.ExitFlag;
      FuEv     = ResultNLP.FuncEv;
      FuncEv   = FuEv;
      GradEv   = ResultNLP.GradEv;
      HessEv   = ResultNLP.HessEv;
      ConstrEv = ResultNLP.ConstrEv;
      % HKH If to run multiMin in case of bad result, delete ~CONVEX
      if any(ExitFlag == [4 10]) & ~CONVEX
         ResultNLP = tomRun('multiMin',ProbNLP,PriLev > 0);
         MM       = ResultNLP.multiMin.Info;
         FuncEv   = FuncEv   + MM.FuncEv;
         GradEv   = GradEv   + MM.GradEv;
         HessEv   = HessEv   + MM.HessEv;
         ConstrEv = ConstrEv + MM.ConstrEv;
      end
   else
      ResultNLP = tomRun('multiMin',ProbNLP,PriLev > 0);
      ExitFlag  = ResultNLP.ExitFlag;
      Inform    = ResultNLP.Inform;
      MM        = ResultNLP.multiMin.Info;
      FuEv      = MM.FuncEv;
      FuncEv    = FuEv;
      GradEv    = MM.GradEv;
      HessEv    = MM.HessEv;
      ConstrEv  = MM.ConstrEv;
      % HKH - Maybe later on try warm start of multiMin with low value of Prob.xInit
      % to assure that global minimum found. In such a case the initial Prob.xInit might have a lower value
      ProbNLP.RandState = [];                      % Avoid using same RandState value again
      % HKH Left this, but could be deleted, because no bug found in SNOPT
      if Inform == 6 & ~strcmpi('npsol',MM.localSolver)
         % Emergency if solver totally fails
         ProbNLP.GO.localSolver = 'npsol';
         ResultNLP = tomRun('multiMin',ProbNLP,PriLev > 0);
         ExitFlag  = ResultNLP.ExitFlag;
         MM        = ResultNLP.multiMin.Info;
         FuEv      = MM.FuncEv;
         FuncEv    = FuncEv   + FuEv;
         GradEv    = GradEv   + MM.GradEv;
         HessEv    = HessEv   + MM.HessEv;
         ConstrEv  = ConstrEv + MM.ConstrEv;
      end
      ProbNLP.xInit     = max(100,10*ProbNLP.N);   % Change default xInit too lower value
   end
   X               = ResultNLP.x_k;
   if isempty(X)
      x            = X;
      Icomp(1)     = 0;
      Ridx         = IntVars;
   else
      [IC,Iidx] = IntComp(X,IntVars,eps_I);
      if length(IC) > 1
         % Select the x with max integer components
         [icMax, idx] = max(IC);
         x            = X(:,idx);
         Ridx         = IntVars(~Iidx(:,idx));
         Icomp(1)     = icMax;
         if idx > 1
            Result.f_k = ResultNLP.multiMin.FX(1,idx);
         end
      else
         x            = X;
         Icomp(1)     = IC;
         Ridx         = IntVars(~Iidx);
      end
   end

   Inform   = ResultNLP.Inform;
   SumMinIt = ResultNLP.Iter;
   v_k      = ResultNLP.v_k;
   if VarSel == 2 & isempty(v_k)
      % Must estimate LM
      % v_k   = nlp_g(x,ProbNLP);
   end
   g_k      = ResultNLP.g_k;
   vIPMin   = []; % Init of Lagrange multipliers in best solution

   SOL=0;
   if CONVEX & isempty(ResultNLP.SOL) 
      if strcmpi(SolverNLP,'minos') | strcmpi(SolverNLP,'sqopt') | ... 
          strcmpi(SolverNLP,'lssol') | strcmpi(SolverNLP,'nlssol') | ...
          strcmpi(SolverNLP,'npsol') | strcmpi(SolverNLP,'snopt') | ...
          strcmpi(SolverNLP,'qpopt') 
         SOL=1;
         ProbNLP.WarmStart=1;
      end
   end

   fNLP0    = ResultNLP.f_k;
   
   if ExitFlag ==4,fNLP0=Inf; end
   
   % Lower bound on the solution
   fLB      = fNLP0;

   if ExitFlag > 0 & ExitFlag ~= 3  %%% & ExitFlag ~= 4
      if PriLev >= 2
         fprintf('No solution found to initial NLP relaxation. ')
         fprintf('ExitFlag %d\n', ExitFlag)
      end
      x(IntVars)          = round(x(IntVars));
      % x(IntVars) = max(x_L(IntVars),min(x_U(IntVars,round(x(IntVars))));
      ResultMIP.x_k       = x;
      ResultMIP.f_k       = fNLP0;
      ResultMIP.v_k       = vIPMin;
      ResultMIP.FuncEv    = FuncEv;
      ResultMIP.GradEv    = GradEv;
      ResultMIP.HessEv    = HessEv;
      ResultMIP.ConstrEv  = ConstrEv;
      ResultMIP.Iter      = 0;
      ResultMIP.ExitFlag  = 4;
      ResultMIP.ExitText  = ExitText(4,WarmStart,1);
      ResultMIP.MinorIt   = SumMinIt;
      ResultMIP.c_k       = [];
      ResultMIP.Ax        = [];
      ResultMIP           = endSolve(Prob,ResultMIP);
      return;
   end

   % Allocate arrays
   if SOL > 0
      XX = spalloc(nb,MxAlloc,min(40000,nb*MaxIter));
      HS = spalloc(nb,MxAlloc,min(40000,nb*MaxIter));
      nS = zeros(1,MxAlloc);
   else
      XX = spalloc(n,MxAlloc,min(40000,n*MaxIter));
      HS = [];
      nS = [];
   end

   % Save solution for 1st node
   if SOL==0
      XX(1:n,1)  = sparse(x);
   else
      XX(1:nb,1) = sparse(ResultNLP.SOL.xs);
      HS(1:nb,1) = sparse(ResultNLP.SOL.hs);
      nS(1,1)    = ResultNLP.SOL.nS;
   end

   % Initialization

   L    = 1;             % Node List
   NODE = 1;

   % Lowest value possible at node, from NLP relaxation
   f_min          = Inf*ones(MxAlloc,1);      
   f_min(1)       = fNLP0;
   % Depth in tree for each node, root node is 0
   Depth          = zeros(MxAlloc,1);      
   % Number of integer components
   Icomp          = zeros(MxAlloc,1);      
   
   pred           = zeros(MxAlloc,1);
   Problem.xVar   = zeros(MxAlloc,1);
   Problem.Bound  = zeros(MxAlloc,1);
   Problem.xBound = zeros(MxAlloc,1);

   SumIter        = 0;
   Iter           = 0;
   fIPIter        = 0;
end                                                   % of COLD START

cpumax = 0;
TIME0  = Prob.TIME0;
BigNum = fLB > 99999 | fLB < -9999;

% Check duality gap
if ~isinf(fIPMin)
    if fIPMin ~= 0
       Gap =  max(0,fIPMin-fLB)/abs(fIPMin);
    else
       Gap =  max(0,fIPMin-fLB);
    end
    if Gap <= DualGap
       GapOK = 1;
    else
       GapOK = 0;
    end
else
    GapOK = 0;
    Gap   = Inf;
end


format compact
while ~(Iter >= MaxIter)
   if isempty(L)
      % There are no more nodes to search !!!
      if fIPMin == Inf % Empty feasible set
         % FAILURE - No feasible solution to MINLP problem
         if PriLev >= 2
            fprintf('\nTotal number of NLP iterations = %d\n',SumMinIt)
            disp('----------------------------------')
            disp('No feasible solution to MINLP problem')
            disp('----------------------------------')
         end
         ResultMIP.x_k       = [];
         ResultMIP.f_k       = fNLP0;
         ResultMIP.v_k       = vIPMin;
         ResultMIP.FuncEv    = FuncEv;
         ResultMIP.GradEv    = GradEv;
         ResultMIP.HessEv    = HessEv;
         ResultMIP.ConstrEv  = ConstrEv;
         Iter                = max(1,Iter);
         ResultMIP.Iter      = Iter;
         ResultMIP.ExitFlag  = 2;
         ResultMIP.ExitText  = ExitText(2,WarmStart,SumIter+Iter);
         ResultMIP.MinorIt   = SumMinIt;
         ResultMIP.DualGap   = Gap;
         ResultMIP.c_k       = [];
         ResultMIP.Ax        = [];
         ResultMIP           = endSolve(Prob,ResultMIP);
         return;
      else
         % BJORN - Typo below?
         % SUCCESS - but NLP failure with have lead to a GAP > 0
         ResultMIP.x_k       = xIPMin;
         ResultMIP.f_k       = fIPMin;
         ResultMIP.v_k       = vIPMin;
         ResultMIP.FuncEv    = FuncEv;
         ResultMIP.GradEv    = GradEv;
         ResultMIP.HessEv    = HessEv;
         ResultMIP.ConstrEv  = ConstrEv;
         Iter                = max(1,Iter);
         ResultMIP.Iter      = Iter;
         ResultMIP.ExitFlag  = 0;
         ResultMIP.ExitText  = ExitText(0,WarmStart,SumIter+Iter);
         ResultMIP.MinorIter = SumMinIt;
         ResultMIP.DualGap   = 0;
         if dLin > 0
            ResultMIP.Ax     = A*xIPMin;
         end
         if dNoL > 0
            %HKH2
            NLP_xc=[]; NLP_c=[]; NARG = [];
            ResultMIP.c_k   = nlp_c(xIPMin,Prob);
         end
         if IterPrint
            PrintIter(2,Iter,[],[],fIPMin,fIPIter,fLB,L,0,0,x,[],[],...
                      Inform,ExitFlag,FuEv,BigNum,'');
            fprintf('--- FilMINT NLP Branch & Bound converged! ')
            fprintf('Its(nodes visited) = %d. ',Iter)
            fprintf('Total NLP Its = %d. ',SumMinIt)
            fprintf('Optimal f(x) =%22.16f.',fIPMin);
            fprintf('\n');
         end
         ResultMIP=endSolve(Prob,ResultMIP);
         return;
      end
   elseif GapOK == 1
      ResultMIP.x_k       = xIPMin;
      ResultMIP.f_k       = fIPMin;
      ResultMIP.v_k       = vIPMin;
      ResultMIP.FuncEv    = FuncEv;
      ResultMIP.GradEv    = GradEv;
      ResultMIP.HessEv    = HessEv;
      ResultMIP.ConstrEv  = ConstrEv;
      Iter                = max(1,Iter);
      ResultMIP.Iter      = Iter;
      ResultMIP.ExitFlag  = 0;
      ResultMIP.ExitText  = ExitText(6,WarmStart,SumIter+Iter);
      ResultMIP.MinorIter = SumMinIt;
      ResultMIP.DualGap   = Gap;
      if dLin > 0
         ResultMIP.Ax     = A*xIPMin;
      end
      if dNoL > 0
         %HKH2
         NLP_xc=[]; NLP_c=[]; NARG = [];
         ResultMIP.c_k   = nlp_c(xIPMin,Prob);
      end
      if IterPrint
         PrintIter(3,Iter,[],[],fIPMin,fIPIter,fLB,L,0,0,x,[],[],...
                   Inform,ExitFlag,FuEv,BigNum,'');
         fprintf('\n--- FilMINT: User defined duality gap reached. ');
         fprintf('Its(nodes visited) = %d. ',Iter)
         fprintf('Total NLP Its =%d. ',SumMinIt)
         fprintf('Optimal f(x) =%22.16f.',fIPMin);
         fprintf('\n');
      end
      ResultMIP          = endSolve(Prob,ResultMIP);
      return;
   end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   Iter = Iter+1;

   % Problem selection and relaxation
   [i,L] = BBnode(NodeSel,L,pred,Depth,Icomp,fIPMin);

   % HKH NOTE - IterPrint > 1 output
   if IterPrint > 1
      xprinti(L,'L:   ',8,14)
      xprinti(Depth(L),'DepL:',8,14)
      xprinti(Icomp(max(1,pred(L)))','Icomp',8,14)
      if BigNum
         xprint(f_min(L),'f_min:',' %8.0f',14)
      else
         xprint(f_min(L),'f_min:',' %8.2f',14)
      end
   end
   Level    = Depth(i);
   
   
   if i == 1 & WarmStart == 0     % 1st step
      x        = full(XX(1:n,1));
      fNode    = fNLP0;
      f_min(1) = fNode;
      
      % Solve CMP, but update LP problem first

      % Make sure derivatives are defined, otherwise evaluate them.
      % evaluate objective derivatives
      size(ProbNLP.A)
      size(ResultNLP.x_k)
      size(ResultNLP.c_k)
      %ax=ProbNLP.A*ResultNLP.x_k
      %ck=ResultNLP.c_k
      %keyboard
      % HKH2 - NOTE!!!!!!!!!!!
      % ResultNLP.x_k may have many columns returned from multiMin
      % Must select ONE. This is done, and saved in XX, and then x = XX above
      % krasch för problem 5, trim loss har många globala optima, fixat bug nedan
      %c_k = [ProbNLP.A*ResultNLP.x_k;ResultNLP.c_k];
      c_k = [ProbNLP.A*x;ResultNLP.c_k];
      g_k = nlp_g(ResultNLP.x_k,Prob);
      % evaluate constraint derivatives
      J_k = [zeros(dLin,n);nlp_dc(ResultNLP.x_k,Prob)];
      ProbLP.x_L = x_L;
      ProbLP.x_U = x_U;
      ProbLP.A(1,1:n) = g_k(:)';
      ProbLP.A(2:m+1,1:n) = J_k;
      ProbLP.b_U(1) = g_k'*ResultNLP.x_k-ResultNLP.f_k;
      ProbLP.b_U(2:m+1) = J_k*ResultNLP.x_k-c_k;
      ResultLP = tomRun(SolverLP, ProbLP, PriLev > 0);
   else
      
      % Reset lower and upper bounds on x
      x_L = x_L0;
      x_U = x_U0;
      j   = i;
      ix  = [];
      % Find the parent node in the tree for the current node, and set the limits in x_L,x_U
      while pred(j) > 1     % Generate NLP relaxation problem
         j  = pred(j);
         ix = [j;ix];
         k  = Problem.xVar(j);
         if Problem.Bound(j)==1
            x_U(k) = min(x_U(k),Problem.xBound(j));
         else
            x_L(k) = max(x_L(k),Problem.xBound(j));
         end
      end
      
      % Solve CMP, but update LP problem first
      
      % Make sure derivatives are defined, otherwise evaluate them.
      % evaluate objective derivatives
      %HKH2
      NARG = []; 
      NLP_x=[]; NLP_f=[]; NLP_xg=[]; NLP_g=[]; 
      NLP_xc=[]; NLP_c=[]; NLP_xdc=[]; NLP_dc=[]; 
      g_k = nlp_g(ResultNLP.x_k,Prob);
      % evaluate constraint derivatives
      J_k = [zeros(dLin,n);nlp_dc(ResultNLP.x_k,Prob)];
      ProbLP.x_L = x_L;
      ProbLP.x_U = x_U;
      ProbLP.A(1,1:n) = g_k(:)';
      ProbLP.A(2:m+1,1:n) = J_k;
      ProbLP.b_U(1) = g_k'*ResultNLP.x_k-ResultNLP.f_k;
      %HKH2 debug printing for bug below
      if 0
      J_k
      m
      n
      xk=ResultNLP.x_k
      ck=ResultNLP.c_k
      keyboard
      end
      %HKH2 BUG! J_k innehåller dLin linjära delar också, men ResultNLP.c_k bara olinjära
      % Krasch för problem 4 och 5
      ProbLP.b_U(2:m+1) = J_k*ResultNLP.x_k-ResultNLP.c_k;
      ResultLP = tomRun(SolverLP, ProbLP, PriLev > 0);

      % Set limits for current node i
      k = Problem.xVar(i);

      if Problem.Bound(i)==1
         x_U(k)=min(x_U(k),Problem.xBound(i));
      else
         x_L(k)=max(x_L(k),Problem.xBound(i));
      end

      % Get optimal x for parent node
      j = pred(i);
      x = full(XX(1:n,j));
      x(k) = max(x_L(k),min(x_U(k),x(k))); % Adjust optimal x for parent node inside x_L,x_U

      if PriLev > 2
         xprint(x,'x:')
      end
  
      if SOL > 0
         ProbNLP.SOL.xs=full(XX(1:nb,j));
         ProbNLP.SOL.hs=full(HS(1:nb,j));
         ProbNLP.SOL.nS=nS(1,j);
         ProbNLP.SOL.xs(k)=x(k);
      else
         ProbNLP.WarmStart=0;
      end

      % A b_L b_U not changed, unless using algorithms squeezing problem using new x_L,x_U
      %ProbNLP.b_L  = b_L; 
      %ProbNLP.b_U  = b_U;
      %ProbNLP.A    = A;
      ProbNLP.x_0  = x;                     % Set x_0 to optimal x for parent node
      ProbNLP.x_L  = x_L; 
      ProbNLP.x_U  = x_U; 
      ProbNLP.QP.DualLimit = fIPMin; % Try stop dual iterations early
      % ProbNLP.P    = [ 'Iter ' num2str(Iter) ' Node ' num2str(i)]; 

      % [ProbNLP.x_L ProbNLP.x_0 ProbNLP.x_U]

      if CONVEX
         % HKH2 - change tomSolve to tomRun, better when both NLP and LP solved
         Result = tomRun(SolverNLP,ProbNLP,PriLev > 0);
         if PriLev > 2
            PrintResult(Result,1);
         end
         FuEv     = Result.FuncEv;
         FuncEv   = FuncEv   + FuEv;
         GradEv   = GradEv   + Result.GradEv;
         HessEv   = HessEv   + Result.HessEv;
         ConstrEv = ConstrEv + Result.ConstrEv;
         if ExitFlag ~= 4  
            % Check if result is infeasible, even if ExitFlag not 4
            % HKH cErrCompute returns || ||_inf
            Result = cErrCompute(Result);
            if Result.bErr > 2*bTol | Result.cErr > 2*cTol, ExitFlag=4; end
         end
         if any(ExitFlag == [4 10]) %%% & ~CONVEX
            % HKH How much work should multiMin do in this case, now just default
            Result   = tomRun('multiMin',ProbNLP,PriLev > 0);
            MM       = Result.multiMin.Info;
            FuEv     = FuEv     + Result.FuncEv;
            FuncEv   = FuncEv   + MM.FuncEv;
            GradEv   = GradEv   + MM.GradEv;
            HessEv   = HessEv   + MM.HessEv;
            ConstrEv = ConstrEv + MM.ConstrEv;
         end
      else
         %HKH Should we utilize new multiMin warm start feature?
         Result   = tomRun('multiMin',ProbNLP,PriLev > 0);
         %Result   = tomRun('multiMin',ProbNLP,2);
         MM       = Result.multiMin.Info;
         FuEv     = MM.FuncEv;
         FuncEv   = FuncEv   + FuEv;
         GradEv   = GradEv   + MM.GradEv;
         HessEv   = HessEv   + MM.HessEv;
         ConstrEv = ConstrEv + MM.ConstrEv;
         ProbNLP.RandState = [];
         if 0 & Result.f_k > fIPMin
            ProbNLP.WarmStart = 1;
            ProbNLP = WarmDefGLOBAL('multiMin',ProbNLP,Result);
            ProbNLP.xInit = 3*Result.multiMin.Info.localTry;
            Result  = tomRun('multiMin',ProbNLP,PriLev > -1);
            ProbNLP.WarmStart = 0;
            MM       = Result.multiMin.Info;
            FuncEv   = FuncEv   + MM.FuncEv;
            GradEv   = GradEv   + MM.GradEv;
            HessEv   = HessEv   + MM.HessEv;
            ConstrEv = ConstrEv + MM.ConstrEv;
         end
      end
      X               = Result.x_k;
      if isempty(X)
         x            = X;
         Icomp(i)     = 0;
         Ridx         = IntVars;
      else
         [IC,Iidx] = IntComp(X,IntVars,eps_I);
         if length(IC) > 1
            % Select the x with max integer components
            [icMax, idx] = max(IC);
            x            = X(:,idx);
            Ridx         = IntVars(~Iidx(:,idx));
            Icomp(i)     = icMax;

            if idx > 1
               Result.f_k = Result.multiMin.FX(1,idx);
            end

         else
            x            = X;
            Icomp(i)     = IC;
            Ridx         = IntVars(~Iidx);
         end
      end
      if IterPrint & ~isnan(Result.f_k)
         %if Result.f_0 > 1E10 | Result.f_k < fLB | Result.f_k > fIPMin
         if Result.f_0 > 1E10 | Result.f_k + abs(Result.f_k)*1E-12 < fLB 
            % xprint(Result.Prob.x_0)
            PrintResult(Result,2);
            if Result.f_k + abs(Result.f_k)*1E-12 < fLB 
               fprintf('\nNOTE - New f_k < fLB \n');
               fprintf('fLB            %30.16f \n', fLB);
               fprintf('Result.f_k     %30.16f \n', Result.f_k);
               fprintf('fLB-Result.f_k %30.16f \n', fLB-Result.f_k);
            end
         end
      end
      v_k      = Result.v_k;
      % HKH VarSel == 2  must be improved
      if VarSel == 2 & isempty(v_k)
         % Must estimate LM
         % v_k   = nlp_g(x,ProbNLP);
      end
      g_k      = Result.g_k;

      SumMinIt = SumMinIt+Result.Iter;
      ExitFlag = Result.ExitFlag;
      Inform   = Result.Inform;

      if ~isempty(x) 
         x   = min(x_U,max(x_L,x));
      end

      if ExitFlag == 0 | ExitFlag == 3
         if SOL == 0 | isempty(Result.SOL)
            XX(1:n,i) = sparse(x);
         else
            XX(:,i)   = sparse(Result.SOL.xs);
            HS(:,i)   = sparse(Result.SOL.hs);
            nS(1,i)   = Result.SOL.nS;
         end

         fNode    = Result.f_k;
         f_min(i) = fNode;
         % Update fLB, lower bound on the solution
         fLB      = min(f_min([i,L]));
      else
         if PriLev > 1
            fprintf('   No Solution, EXIT = %4.0f\n',ExitFlag)
         end
         fNode    = Inf;
         f_min(i) = Inf;
      end
   end
   ENDTREE = 1;
   xC      = [];
   xI      = [];
   DelL    = [];
   if fNode <  fIPMin
      ENDTREE=0;
      % Check if solution is integer
      if nI == IC
         IntVar=1;
      else
         % Branch on some non-integer variable
         if VarSel == 1
            [IntVar,xC,xI] = BranchVar(VarSel,x,Ridx,VarWeight,eps_I,PriLev);
         elseif VarSel == 2
            if isempty(v_k), v_k = VarWeight; end
            if isempty(v_k), v_k = rand(n,1); end % Panic right now, must estimate v_k
            [IntVar,xC,xI] = BranchVar(VarSel,x,Ridx,v_k(1:n),eps_I,PriLev);
         end
      end
      if IntVar
         if fNode <  fIPMin  % Check if new integer value lower than best found
            fIPMin           = fNode;
            fIPIter          = Iter;
            x(IntVars)       = round(x(IntVars));
            xIPMin           = x(1:n);
            vIPMin           = v_k;
            DelL             = f_min(L) >= fIPMin;
            XX(:,L(DelL))    = 0;
            if SOL > 0
               HS(:,L(DelL)) = 0;
               nS(1,L(DelL)) = 0;
            end
            L                = L(~DelL);

            if IterPrint
               PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),...
                         nI-Icomp(i),x,xC,xI, ...
                         Inform,ExitFlag,FuEv,BigNum,'End of tree. Found IP: ');
            end
            if PriLev > 2
               xprint(xIPMin,'x:');
            end

            if abs(fIPMin - nlp_f(xIPMin,ProbNLP)) < 1E-10 %CHANGED
               ENDTREE=1;
            else
               xC=[];
            end
            if fIPMin ~= 0
               Gap =  max(0,fIPMin-fLB)/abs(fIPMin);
            else
               Gap =  max(0,fIPMin-fLB);
            end
            if Gap <= DualGap
               ENDTREE=1;
               GapOK = 1;
            end
         end
      else
         % Apply heuristics
         RndHeu = 0;
         if ROUNDH                 % Check rounding heuristic
            xRS          = x;
            xRS(IntVars) = round(x(IntVars));
            if dLin > 0
               % HKH consViol returns || ||_1, max(LVi) = || ||_inf
               % HKH consViol more efficient than cErrCompute
               [LV LVi]     = consViol(ProbNLP.A*xRS,b_L,b_U,absviol);
            else
               LV           = 0;
            end
            if LV <= bTol
               fRS                    = nlp_f(xRS,ProbNLP);
               if fRS < fIPMin
                  if dNoL > 0
                     %HKH2
                     NLP_xc=[]; NLP_c=[]; NARG = [];
                     c_k              = nlp_c(xRS,ProbNLP);
                     [NV NVi]         = consViol(c_k,c_L,c_U,absviol);
                  else
                     NV               = 0;
                  end
                  if NV <= dNoL*cTol
                     RndHeu           = 1;
                     fIPMin           = fRS;
                     fIPIter          = Iter;
                     xIPMin           = xRS;
                     vIPMin           = v_k;
                     DelL             = f_min(L) >= fIPMin;
                     XX(:,L(DelL))    = 0;
                     if SOL > 0
                        HS(:,L(DelL)) = 0;
                        nS(1,L(DelL)) = 0;
                     end
                     L                = L(~DelL);

                     if IterPrint
                        PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,...
                                  sum(DelL),nI-Icomp(i),x,xC,xI,Inform,ExitFlag,...
                                  FuEv,BigNum,'Rounding Heuristic: ');
                     end
                  end
               end
            end
         end
         if KNAPSACK > 0 & RndHeu == 0
            % Apply rounding knapsack heuristic
            xKS          = x;
            xKS(IntVars) = floor(x(IntVars)+eps_I*max(1,abs(x(IntVars))));
            if dLin > 0
               if nMK > 0
                  xKS(xMK) = 0;                               % Set slacks to 0
                  xKS(xMK) = b_U(bMK) - ProbNLP.A(bMK,:)*xKS; % Find slack values
                  xKS(xIK) = round(xKS(xIK));                 % Must round integer slacks
                  xKS(xMK) = max(x_L(xMK),min(x_U(xMK),xKS(xMK))); % Inside bounds
               end
               [LV LVi]     = consViol(ProbNLP.A*xKS,b_L,b_U,absviol);
            else
               LV           = 0;
            end
            if LV < bTol
               fKS                    = nlp_f(xKS,ProbNLP);
               if fKS < fIPMin
                  % Also check nonlinear constraints
                  if dNoL > 0
                     %HKH2
                     NLP_xc=[]; NLP_c=[]; NARG = [];
                     c_k              = nlp_c(xKS,ProbNLP);
                     [NV NVi]         = consViol(c_k,c_L,c_U,absviol);
                  else
                     NV               = 0;
                  end
                  if NV <= dNoL*cTol
                     % Accept new point
                     fIPMin           = fKS;
                     fIPIter          = Iter;
                     xIPMin           = xKS;
                     vIPMin           = v_k;
                     DelL             = f_min(L) >= fIPMin;
                     XX(:,L(DelL))    = 0;
                     if SOL > 0
                        HS(:,L(DelL)) = 0;
                        nS(1,L(DelL)) = 0;
                     end
                     L                = L(~DelL);
                     if IterPrint
                        PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),...
                                  nI-Icomp(i),x,xC,xI,...
                                  Inform,ExitFlag,FuEv,BigNum,'Knapsack Heuristic: ');
                     end
                  end
               end
            end
         end
      end
   end
   if fNode > fIPCut              % Cut the tree, user defined fIPCut
      ENDTREE=1;
   end
   if ~ENDTREE & ~isempty(xC)     % Make divisions
      NODE                        = NODE + 1;
      Problem.xVar(NODE:NODE+1)   = xC;
      Problem.Bound(NODE:NODE+1)  = [-1; 1];
      Problem.xBound(NODE:NODE+1) = [xI+1; xI];
      L                           = [L NODE NODE+1];
      f_min(NODE:NODE+1)          = fNode;
      pred(NODE:NODE+1)           = i;
      Depth(NODE)                 = Level + 1;
      NODE                        = NODE  + 1;
      Depth(NODE)                 = Level + 1;
      if PriLev > 3
         Ntree    = [i;NODE];
         j        = i;
         while pred(j) > 0        % Generate tree
            j     = pred(j);
            Ntree = [j;Ntree];
         end
         disp('  # Node Pred xVar xBound Bound           Value')
         fprintf('%3d %4d %4d %4d\n',1,Ntree(1), pred(1),Problem.xVar(1));
         for ii = 2:length(Ntree)
             jj = Ntree(ii);
             kk = Problem.xVar(jj);
             fprintf('%3d %4d %4d %4d %6d %5d %15.7f\n',ii,jj,...
                  pred(jj),kk,Problem.xBound(jj),Problem.Bound(jj),x(kk));
         end
      end

      if PriLev > 4
         % Dump of all nodes
         fprintf('  Node  Pred  xVar xBound Bound Depth IntComp      f_min\n')
         fprintf('%5d %5d %5d  %5d %5d %5d %5d %15.8f\n',[[1:NODE]' pred(1:NODE) ...
                 Problem.xVar(1:NODE) Problem.xBound(1:NODE)...
                 Problem.Bound(1:NODE) Depth(1:NODE) Icomp(1:NODE) f_min(1:NODE) ]');
      end
   end
   if IterPrint
      PrintIter(1,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),nI-Icomp(i),x,xC,xI,...
                Inform,ExitFlag,FuEv,BigNum,'');
   end

   if mod(Iter,SaveFreq) == 0
      time = fix(clock);
      try
         TotIter   = SumIter+Iter;
         RandState = rand('state');
         save('FilMINTSave.mat','Name','xIPMin','fIPMin','fIPIter','vIPMin',...
                            'FuncEv','GradEv','HessEv', 'ConstrEv', ...
                            'XX','HS','nS', ...
                            'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                            'Depth','Icomp',...
                            'Gap','SumMinIt','RandState','time')
      catch
         warning('Failed to save warm start information to FilMINTSave.mat');
         err = lasterror;
         disp(err.message);
         disp('Warm start information is available in Result.MINLP');
      end
   end
   if wait 
         disp('Press any key to continue ...')
         pause(1)
   end
end % while

ResultMIP.x_k                  = xIPMin;
ResultMIP.f_k                  = fIPMin;
ResultMIP.v_k                  = vIPMin;
ResultMIP.FuncEv               = FuncEv;
ResultMIP.GradEv               = GradEv;
ResultMIP.HessEv               = HessEv;
ResultMIP.ConstrEv             = ConstrEv;

% Warm start info saved
ResultMIP.MINLP.Name      = Name;
ResultMIP.MINLP.xIPMin    = xIPMin;
ResultMIP.MINLP.fIPMin    = fIPMin;
ResultMIP.MINLP.fIPIter   = fIPIter;
ResultMIP.MINLP.vIPMin    = vIPMin;
ResultMIP.MINLP.FuncEv    = FuncEv;
ResultMIP.MINLP.GradEv    = GradEv;
ResultMIP.MINLP.HessEv    = HessEv;
ResultMIP.MINLP.ConstrEv  = ConstrEv;
ResultMIP.MINLP.XX        = XX(:,1:NODE);
if SOL > 0
   ResultMIP.MINLP.HS     = HS(:,1:NODE);
   ResultMIP.MINLP.nS     = nS(1,1:NODE);
else
   ResultMIP.MINLP.HS     = HS;
   ResultMIP.MINLP.nS     = nS;
end
ResultMIP.MINLP.L         = L;
ResultMIP.MINLP.NODE      = NODE;

TotIter                        = SumIter+Iter;

ResultMIP.MINLP.TotIter   = TotIter;
ResultMIP.MINLP.f_min     = f_min(1:NODE);
ResultMIP.MINLP.pred      = pred(1:NODE);
ResultMIP.MINLP.Problem.xVar     = Problem.xVar(1:NODE);
ResultMIP.MINLP.Problem.Bound    = Problem.Bound(1:NODE);
ResultMIP.MINLP.Problem.xBound   = Problem.xBound(1:NODE);
% ResultMIP.MINLP.Problem = Problem;
ResultMIP.MINLP.Depth     = Depth;
ResultMIP.MINLP.Icomp     = Icomp;
ResultMIP.MINLP.Gap       = Gap;
ResultMIP.MINLP.SumMinIt  = SumMinIt;
RandState                      = rand('state');
ResultMIP.MINLP.RandState = RandState;
time                           = fix(clock);
ResultMIP.MINLP.time      = time;

try
   save('FilMINTSave.mat','Name', 'xIPMin','fIPMin','fIPIter','vIPMin',...
                      'FuncEv','GradEv','HessEv', 'ConstrEv', ...
                      'XX','HS','nS', ...
                      'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                      'Depth','Icomp',...
                      'Gap','SumMinIt','RandState','time');
catch
   warning('Failed to save warm start information to FilMINTSave.mat');
   err = lasterror;
   disp(err.message);
   disp('Warm start information is available in Result.MINLP');
end

if cpumax
   ResultMIP.ExitFlag=99;
   ResultMIP.ExitText=ExitText(9,WarmStart,TotIter);
   if PriLev >= 1
      fprintf('\n--- MAXCPU exceeded. Branch and Bound iterations ITER = %d\n',Iter)
      fprintf('\n--- Total number of NLP iterations = %d\n',SumMinIt)
   end
else
   ResultMIP.ExitFlag=1;
   ResultMIP.ExitText=ExitText(1,WarmStart,TotIter);
   if PriLev >= 1
      fprintf('\n--- TOO MANY Branch and Bound ITERATIONS. ITER = %d\n',Iter)
      fprintf('\n--- Total number of NLP iterations = %d\n',SumMinIt)
   end
end
ResultMIP.Iter      = Iter;
ResultMIP.MinorIter = SumMinIt;
ResultMIP.DualGap   = Gap;
ResultMIP           = endSolve(Prob,ResultMIP);

% ------------------------------
function Result = cErrCompute(Result)
% ------------------------------
x_k  = Result.x_k;
if isempty(x_k)
   Result.bErr = Inf;
   Result.cErr = Inf;
   return
end
Prob = Result.Prob;
bErr = 0;
cErr = 0;
if Prob.mLin > 0
   if isempty(Result.Ax)
      Ax = Result.A * x_k;
      Result.Ax = Ax;
   else
      Ax = Result.Ax;
   end
%[Prob.b_L Ax Prob.b_U]
   ixE  = find(Prob.b_L == Prob.b_U);
   z =max(1,abs(Prob.b_U));
   if ~isempty(ixE)
      bErr = max(bErr, max(abs(Prob.b_U(ixE)-Ax(ixE))./z(ixE)));
   end
   ixU  = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_U));
   if ~isempty(ixU)
      bErr = max(bErr, max(max(0,Ax(ixU)-Prob.b_U(ixU))./z(ixU)));
   end
   ixL  = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_L));
   if ~isempty(ixL)
      z =max(1,abs(Prob.b_L(ixL)));
      bErr = max(bErr, max(max(0,Prob.b_L(ixL)-Ax(ixL))./z));
   end
end
if Prob.mNonLin > 0
   c_k  = Result.c_k;
%[Prob.c_L c_k Prob.c_U]
   ixE  = find(Prob.c_L == Prob.c_U);
   z =max(1,abs(Prob.c_U));
   if ~isempty(ixE)
      cErr = max(cErr, max(abs(Prob.c_U(ixE)-c_k(ixE))./z(ixE)));
   end
   ixU  = find(Prob.c_L ~= Prob.c_U & ~isinf(Prob.c_U));
   if ~isempty(ixU)
      cErr = max(cErr, max(max(0,c_k(ixU)-Prob.c_U(ixU))./z(ixU)));
   end
   ixL  = find(Prob.c_L ~= Prob.c_U & ~isinf(Prob.c_L));
   if ~isempty(ixL)
      z =max(1,abs(Prob.c_L(ixL)));
      cErr = max(cErr, max(max(0,Prob.c_L(ixL)-c_k(ixL))./z));
   end
end
Result.bErr = bErr;
Result.cErr = cErr;


% -----------------------------------------------
function [h_L1, h] = consViol(z,L,U,absviol)
% -----------------------------------------------
if absviol == 1
   h =-min(0,z-L)-min(0,U-z);
else
   h =-min(0,(z-L)./max(1,abs(L)))-min(0,(U-z)./max(1,abs(U)));
end
h_L1 = sum(h);


% ------------------------------
function Text = ExitText(Inform,WarmStart,Its)
% ------------------------------

switch  Inform
   case 0
     if WarmStart
        Text = ['Optimal solution found. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'Optimal solution found';
     end
   case 1
     if WarmStart
        Text = ['Max iterations. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'Maximal number of iterations reached';
     end
   case 2
     if WarmStart
        Text = ['Empty feasible set, no integer solution found. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'Empty feasible set, no integer solution found';
     end
   case 4
     if WarmStart
        Text = ['No solution found to NLP relaxation. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'No solution found to NLP relaxation';
     end
   case 5
     if WarmStart
        Text = ['Illegal x_0 found in NLP relaxation. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'Illegal x_0 found in NLP relaxation';
     end
   case 6
     if WarmStart
        Text = ['User defined duality gap reached. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'User defined duality gap reached';
     end
   case 9
     if WarmStart
        Text = ['User defined maximal CPU Time used. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'User defined maximal CPU Time used';
     end
end

% -----------------------------------------------
function PrintIter(Type,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,DelL,nInts,x,xC,xI,...
                   Inform,ExitFlag,FuEv,BigNum,Txt1);
% -----------------------------------------------
if Type < 4
   switch Type
   case 1
     fprintf('It%6.0f:',Iter);
     fprintf('%3d ',Level);
   case 2
      fprintf('It%6.0f:L[] ',Iter);
   case 3
      fprintf('It%6.0f:Gap ',Iter);
   end

   if Type == 1
      if BigNum
         fprintf('fNLP %9.0f ',fNode);
      else
         fprintf('fNLP %9.3f ',fNode);
      end
   else
      fprintf('%s',blanks(15));
   end
else
   fprintf('%s%s',blanks(28-length(Txt1)),Txt1);
end

fprintf('fIP %10.3f ',fIPMin);
if isinf(fIPMin)
   fprintf('%s',blanks(7));
else
   fprintf('at%4d ',fIPIter);
end

if Type < 4
   if BigNum
      fprintf('LowB%10.0f ',fLB);
   else
      fprintf('LowB%10.3f ',fLB);
   end
   fprintf('Gap%8.2f ',max(0,fIPMin-fLB));
   if fIPMin ~= 0 & ~isinf(fIPMin)
      fprintf('%8.3f%% ',100*max(0,fIPMin-fLB)/fIPMin);
   else
      fprintf('%s',blanks(10));
   end
else
   fprintf('%s',blanks(37));
end
if Type == 1 | Type == 4
   fprintf('|L|%4d ',length(L));
   if isempty(x)
      fprintf('%s',blanks(5));
   else
      fprintf('I%3d ',nInts);
   end
   if DelL == 0
      fprintf('%s',blanks(5));
   else
      fprintf('D%3d ',DelL);
   end
end
if Type == 1
   if isempty(xC)
      fprintf('%s',blanks(11));
   else
      fprintf('xC%3d ',xC);
      fprintf('xI%2d ',xI);
   end
   fprintf('! Info%2d ',Inform);
   fprintf('Ex%2d ',ExitFlag);
   fprintf('FuEv%6d',FuEv);

   time      = fix(clock);
   tt = time([3 2 4 5]);
   if tt(4) < 10
      fprintf(' %d/%d %d:0%1d ', tt);
   else
      fprintf(' %d/%d %d:%2d ', tt);
   end
end
fprintf('\n');

% MODIFICATION LOG
%
% 091118  hkh  Written, based on minlpSolve
