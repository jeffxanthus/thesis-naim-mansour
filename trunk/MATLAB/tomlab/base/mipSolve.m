% mipSolve.m:
%
% function Result = mipSolve(Prob)
%
% Branch & Bound algorithm for Mixed-Integer Programming (MIP)
% using LP relaxation (Formulated as min-IP)
%
% Solving MIP on the LP form with generalized lower and upper bounds:
%
%        min    c' * x.  x in R^n, b_L <= A x <= b_U, x_L <= x <= x_U
%
% Any of the x could be set as integer valued
%
% Note that the TOMLAB routine cpTransf.m can be used to rewrite
% an LP problem on another form to this standard form
%
% INPUT PARAMETERS
% Fields in Prob:
%   c:      The vector c in c'x.
%   A:      The linear constraint matrix
%   b_L:    Lower bounds on linear constraints.
%           If empty, assumed to be == b_U
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x. If empty assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%
%   x_0:    Starting point x (If EMPTY, the LP solver solves a Phase I LP
%           to find a feasible point. Some LP solvers chooses own x_0).
% MaxCPU:   Maximal CPU Time (in seconds) to be used by mipSolve,
%           stops with best point found
%
% PriLev    Print level in mipSolve (default 1). Also see optParam.IterPrint.
% PriLevOpt Printing level in subsolvers, see the documentation for the solver
%
% WarmStart If true, >0, mipSolve reads the output from the last run from Prob.mipSolve, 
%           if it exists.  If it doesn't exist, mipSolve attempts to open and read warm 
%           start data from mat-file mipSolveSave.mat. 
%           mipSolve uses the warm start information to continue from the last run.
%           The mat-file mipSolveSave.mat is saved every Prob.MIP.SaveFreq iteration
%
%
% QP.B:     Active set B_0 at start. (Only used by lpSimple and DualSolve)
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound
%           -1 = Variable x(i) is set on its upper bound. NOT USED
%           If EMPTY, lpSimplex finds active set (if used for LP solution).
%
% SolverLP  Name of the solver used for the initial LP subproblem. If empty,
%           picked from a list, best available with a license
%           Note - If using MINOS, use the special LP interface LP-MINOS
% SolverDLP Name of the solver used for dual LP subproblems. If empty, SolverLP is used.
%           If SolverLP/SolverDLP is a SOL solver (SNOPT, MINOS or NPSOL), 
%           the SOL.optPar and SOL.PrintFile is used:
%           See help minosTL.m, npsolTL.m or snoptTL.m for how to set these parameters
%
% ---------------------------------------
% MIP       Structure in Prob, Prob.MIP
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
%           Prob.MIP.VarWeight = c; for knapsack problems improve convergence.
% DualGap   mipSolve stops if the duality gap is less than DualGap
%           E.g. DualGap = 1, implies a stop at first integer solution
%           E.g. DualGap = 0.01, stop if solution < 1% from optimal solution
% fIP       An upper bound on the IP value wanted. Makes it possible to
%           cut branches and avoid node computations.  Used even if xIP not given
% xIP       The x-values giving the fIP value, if a solution (xIP,fIP) is known.
%
% NodeSel   Node selection method in branch and bound
%           = 0 Depth First. Priority on  nodes with more integer components.
%           = 1 Breadth First. Priority on  nodes with more integer components.
%           = 2 Depth First. When integer solution found, use NodeSel = 1 (default)
%           = 3 Pure LIFO (Last in, first out) Depth First
%           = 4 Pure FIFO (First in, first out) Breadth First
%           = 5 Pure LIFO Depth First. When integer solution found, use NodeSel 4
%
% VarSel    Variable selection method in branch and bound
%           = 1 Use variable with most fractional value
%           = 2 Use gradient and distance to nearest integer value 
% KNAPSACK  If = 1, use a knapsack heuristic. Default 1.
% ROUNDH    If = 1, use a rounding heuristic. Default 1.
%
% SaveFreq  Warm start info saved on mipSolveSave.mat every SaveFreq iteration.
%           (default -1, i.e. no warm start info is saved)
% ----------------------------------------------------------------------
%
% Options for the LP subproblem
% Prob.Solver.Method  Rule to select new variables in DualSolve/lpSimplex:
%           Note - this option is not used for other LP solvers.
%           = 0 Minimum Reduced Cost, sort variables increasing. (DEFAULT)
%           = 1 Bland's Anti-cycling Rule
%           = 2 Minimum Reduced Cost, Dantzig's rule
%
% -------------------------------------------------------------------------------
% optParam  Structure in Prob. Fields used in Prob.optParam, also in sub solvers:
% -------------------------------------------------------------------------------
% MaxIter   Maximal number of iterations. max(10*dim(x),100) is DEFAULT.
% IterPrint Print short information each iteration (PriLev > 0 ==> IterPrint = 1)
%           One line each iteration with: Iteration number: Depth in tree (symbol 
%           L[] - empty list, symbol Gap - Dual Gap convergence; date/time stamp, 
%           fLP (Optimal f(x) current node), fIPMin (Best integer feasible f(x) found)
%           LowBnd (Lower bound on optimal integer feasible f(x)),
%           Dual Gap in absolut value and percent.
%            The length of the node list L, |L|
%           The Inform and ExitFlag the solver returned at the current node
% bTol      Linear constraint violation convergence tolerance
%
% Fields in Prob.optParam that are used in lpSimplex and DualSolve (in case they are used)
% eps_f     Tolerance used for function value tests
% eps_Rank  Rank tolerance
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
%            == 4  => No feasible point found running LP relaxation
%            == 5  => Illegal x_0 found in LP relaxation
%            == 99 => Maximal CPU Time used (cputime > Prob.MaxCPU)
%   ExitTest Text string giving ExitFlag and Inform information
%   Inform   If ExitFlag > 0, Inform=ExitFlag, except if
%            duality gap less than given tolerance (Inform = 6)
%   DualGap  Relative duality gap, max(0,fIPMin-fLB)/|fIPMin|, if fIPMin ~=0;
%            max(0,fIPMin-fLB) if fIPMin == 0. If fIPMin ~=0:
%            Scale with 100, 100*DualGap, to get the percentage duality gap.
%            For absolute value duality gap: scale with fIPMin, fIPMin * DualGap
%   x_k      Solution
%   v_k      Lagrange parameters (n+m) vector. [Bounds (n); Linear constraints (m)]
%   QP.B     B  Optimal set. B(i)==0, include variable x(i) in basic set.
%            sum(B==1)==length(b)  holds. See QP.B as input.
%   QP.y     Dual parameters y (also part of v_k)
%   f_k      Function value c'*x, fIPMin (Best integer feasible solution)
%   g_k      Gradient c
%   Solver   mipSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%
%   mipSolve 
%       A structure with warm start information. Use with WarmDefGLOBAL, see example below
%
% -- Calling mipSolve
%
% Recommended is to first set IterPrint, to get information each iteration
%      Prob.optParam.IterPrint = 1;
%
% Driver call, including printing with level 2:
%      Result = tomRun('mipSolve',Prob,2);
%
% Direct solver call:
%      Result = mipSolve(Prob);
%      PrintResult(Result);
%
% -- Warm start
% To make a restart (warm start), just set the warm start flag, and call
% mipSolve once again:
%
%      Prob.WarmStart = 1;
%      Result         = tomRun('mipSolve', Prob, 2);
%
% mipSolve will read warm start information from the mipSolveSave.mat file
%
% Another warm start (with same MaxFunc) is made by just calling tomRun again
%      Result = tomRun('mipSolve', Prob, 2);
%
% To make a restart from the warm start information in the Result
% structure, make a call to WarmDefGLOBAL before calling mipSolve.
% WarmDefGLOBAL moves information from the Result structure
% to the Prob structure and sets the warm start flag, Prob.WarmStart = 1; 
%
%      Prob = WarmDefGLOBAL('mipSolve', Prob, Result);
% where Result is the result structure returned by the previous run.
% A warm start (with same MaxIter) is done by just calling tomRun again
%      Result = tomRun('mipSolve', Prob, 2);
%
% To make another warm start with new MaxIter 100, say, redefine MaxIter as:
%      Prob.optParam.MaxIter = 100;
% Then repeat the two lines:
%      Prob = WarmDefGLOBAL('mipSolve', Prob, Result);
%      Result = tomRun('mipSolve', Prob, 2);
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1994-2010 by Tomlab Optimization Inc., $Release: 7.5.0$
% Written Nov 30, 1994.   Last modified June 5, 2010.

function ResultMIP = mipSolve(Prob)

if nargin < 1
   error('mipSolve needs input structure Prob');
end

solvType  = checkType('mip');
Prob      = ProbCheck(Prob,'mipSolve',solvType);
Prob      = iniSolveMini(Prob);
Prob.QP.F = [];
absviol   = 0;   % Check relative error in linear constraints

MaxCPU    = DefPar(Prob,'MaxCPU',Inf);
dLin      = size(Prob.A,1);
n         = Prob.N;
nb        = n+dLin+1;
x         = Prob.x_0(:);  

if any(isinf(x) | isnan(x)), x=[]; end
if isempty(x), x = zeros(n,1); end

x_L       = Prob.x_L(:);
if isempty(x_L), x_L=zeros(n,1); end
x_U       = Prob.x_U(:);
if isempty(x_U), x_U=Inf*ones(n,1); end
% Safe-guard starting point
x         = max(x_L(1:n),min(x,x_U(1:n)));

Name      = deblank(Prob.Name);
B         = Prob.QP.B(:);
x_L0      = x_L;
x_U0      = x_U;
b_L       = Prob.b_L;
b_U       = Prob.b_U;
Prob.x_L  = x_L;
Prob.x_U  = x_U;
optParam  = Prob.optParam;
MaxIter   = optParam.MaxIter;     % Maximal number of iterations
eps_f     = optParam.eps_f;    
bTol      = optParam.bTol;        % Linear constraint violation convergence tolerance
MxAlloc   = MaxIter*2+2;

PriLev    = Prob.PriLev;          % Print level
PriLevOpt = Prob.PriLevOpt;       % Print level in LP and Dual LP solver
IterPrint = optParam.IterPrint;   % Print short information each iteration

if isempty(PriLev), PriLev    = 0; end
if PriLev > 0,      IterPrint = 1; end
 
% Integer variables
IntVars   = DefPar(Prob.MIP,'IntVars',[]);
% Logical vector for integers
IV        = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('mipSolve: Illegal IntVars vector');
   end
   IV(IntVars)           = 1;
end
IntVars          = find(IV);
Prob.MIP.IntVars = IntVars;
nI               = length(IntVars);
eps_I            = 1E-8;             % Tolerance when x_i is considered integer valued

ResultMIP = ResultDef(Prob);

fIP       = DefPar(Prob.MIP,'fIP',[]);
fIPMin    = fIP;
xIPMin    = DefPar(Prob.MIP,'xIP',[]);

if isempty(fIPMin)
   fIPMin = Inf;                   % Current BEST integer solution
   if ~isempty(xIPMin)             % User has supplied an x
      if length(xIPMin) < n
         xIPMin       = [xIPMin;zeros(n-length(xIPMin))];
      end
      xIPMin(IntVars) = round(xIPMin(IntVars));
      fIPMin          = c'*xIPMin;  % Compute the function value at the given point
   end
else
   if isempty(xIPMin)
      fIPMin = Inf;                 % Current BEST integer solution
   else
      xIPMin(IntVars) = round(xIPMin(IntVars));
   end
end

if ~isinf(fIPMin) 
   if Prob.mLin > 0
      [LV LVi]     = consViol(Prob.A*xIPMin,b_L,b_U,absviol);
   else
      LVi          = 0;
      LV           = 0;
   end
   XLV             = max(0,x_L-xIPMin);
   XUV             = max(0,xIPMin-x_U);
   if any(XLV > 1E-10  | XUV > 1E-10) | any(LVi > bTol) 
      if IterPrint | PriLev > 0
         fprintf('Sum of error in linear constraints     %f\n',LV);
         fprintf('Sum of error in lower bound violations %e\n',sum(XLV));
         fprintf('Sum of error in upper bound violations %e\n',sum(XUV));
         fprintf('xIPMin is not feasible\n');
      end
      fIPMin       = Inf;     
      xIPMin       = [];     
   end
end

if ~isempty(fIPMin) & isempty(xIPMin)
   fIPCut = fIP;
else
   fIPCut = Inf;
end

if isfield(Prob.MIP, 'DualGap')
   if isempty(Prob.MIP.DualGap)
      DualGap = 0;
   else
      DualGap = Prob.MIP.DualGap;
   end
else
   DualGap = 0;
end
BIPMin = [];
yIPMin = [];
vIPMin = [];

if ~isfield(Prob.MIP,'VarWeight')
   Prob.MIP.VarWeight = [];
   VarWeight = [];
else
   VarWeight = Prob.MIP.VarWeight(:);
end

% Warm start info saving frequency
SaveFreq           = DefPar(Prob.MIP,'SaveFreq',-1);

% Node selection method
NodeSel            = DefPar(Prob.MIP,'NodeSel',2);

% Variable selection method
VarSel            = DefPar(Prob.MIP,'VarSel',1);

% Heuristics
KNAPSACK          = DefPar(Prob.MIP,'KNAPSACK',1);
ROUNDH            = DefPar(Prob.MIP,'ROUNDH',1);

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

ResultMIP                 = ResultDef(Prob);
ResultMIP.Solver          = 'mipSolve';
ResultMIP.Prob            = Prob;
ResultMIP.mipSolve        = [];

% Index selection rule in DualSolve and lpSimplex
Method             = Prob.Solver.Method;
if isempty(Method),  Method = 0; end
Prob.Solver.Method = Method;

if isempty(Prob.SolverLP)
   if checkLicense('minos')
      SolverLP='lp-minos';
   elseif checkLicense('bqpd')
      SolverLP='bqpd';
   elseif checkLicense('xa')
      SolverLP='xa';
   else
      SolverLP='milpSolve';
   end
   % SolverLP=GetSolver('lp',Prob.LargeScale);
   Prob.SolverLP = SolverLP;
else
   SolverLP=Prob.SolverLP;
end

ResultMIP.SolverAlgorithm = 'Branch & Bound MILP - ';

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
case 6
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' Depth First, then improve bounds.'];
otherwise
   NodeSel = 2;
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' Depth First, then Breadth.'];
end

if ~isempty(VarWeight)
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' Priority weights.'];
end

if PriLev > 0
   fprintf('%s',ResultMIP.SolverAlgorithm);
   fprintf('\n');
end

if isempty(Prob.SolverDLP)
   SolverDLP      = SolverLP;
   Prob.SolverDLP = SolverDLP;
else
   SolverDLP      = Prob.SolverDLP;
end

ResultMIP.f_0  = [];
ResultMIP.x_k  = [];
ResultMIP.f_k  = [];
ResultMIP.v_k  = [];
ResultMIP.QP.B = [];
ResultMIP.QP.y = [];

% Setup Structure used in 1st LP call
ProbLP      = CreateProbQP(Prob, 'lp', max(5000,10*n), PriLev-3, 0);
ProbLP.QP.c = Prob.QP.c;
ProbLP.x_L  = x_L; 
ProbLP.x_U  = x_U; 
ProbLP.P    = 0; 
ProbLP.SOL  = Prob.SOL;
ProbLP.x_0  = x;

if PriLev ==2 & PriLevOpt >=2
   ProbLP.PriLevOpt   = PriLevOpt-1;
else
   ProbLP.PriLevOpt   = max(0,PriLevOpt-1);
end

% Setup Structure used in dual LP calls
ProbDLP = ProbLP;

% mipSolve warm start parameters
MILP = DefPar(Prob,'mipSolve',[]);

if Prob.WarmStart >= 1                    % Restart with values from previous run.

   Name1 = Name;  % Problem name
   
   % If the WarmStart structure exists, then use that warm start data.
   if ~isempty(MILP)
      Name = MILP.Name;
      if strcmp(Name1,Name)
         if PriLev >= 2
            fprintf(['Restart with %d iterations from last ' ...
                  'runs.\nUsing warm start data from structure: ' ...
                  'Prob.mipSolve.\n'],...
                  MILP.TotIter);
         end
      else
         Prob.WarmStart = 0;
         if PriLev >= -1000
            fprintf('Previous run was with Problem %s\n',Name);
            fprintf('This run is with Problem %s\n',Name1);
            fprintf(['Impossible to do restart using data from ' ...
                     'structure: Prob.mipSolve.\n']);
         end
      end

   else % If no WarmStart structure exists, look for a file.
      filename = 'mipSolveSave.mat';
     
      if exist(filename,'file')~=2
         fprintf(['Couldn''t find warm start info structure nor warm ' ...
                  'start file in path.\nNo warm start possible.\n']);
         Prob.WarmStart = 0;
      else
         load(filename,'Name');

         if strcmp(Name1,Name)
            MILP = load(filename,'xIPMin','fIPMin','fIPIter','vIPMin',...
                            'XX','HS','BB', ...
                            'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                            'Depth','Icomp',...
                            'Gap','SumMinIt','fDual','time');
        
            if PriLev >= 2
               fprintf(['Restart with %d iterations from last ' ...
                       'runs.\nUsing warm start data from file: %s.\n'],MILP.TotIter, filename);
            end
         
         else
            Prob.WarmStart = 0;
            if PriLev >= -1000
               fprintf('Previous run was with Problem %s\n',Name);
               fprintf('This run is with Problem %s\n',Name1);
               fprintf(['Impossible to do restart using data from file: %s.\n'], filename);
            end
        end
     end
   end
end

WarmStart           = Prob.WarmStart;

if WarmStart == 2 
   % Shrinking the number of nodes too costly now
   % Get information from last run
   xIPMin    = MILP.xIPMin;
   fIPMin    = MILP.fIPMin;
   fIPIter   = MILP.fIPIter;
   vIPMin    = MILP.vIPMin;
   L         = MILP.L;
   NODE      = MILP.NODE;
   TotIter   = MILP.TotIter;
   Gap       = MILP.Gap; % Not needed, recomputed below
   SumMinIt  = MILP.SumMinIt;
   fDual     = MILP.fDual;
   time      = MILP.time ;
   if IterPrint > 0
      fprintf('Warm Start with %d nodes. ',NODE);
   end
   if isempty(MILP.HS)
      SOL               = 0;
   else
      SOL               = 1;
      ProbNLP.WarmStart = 1;
   end
   % k is the number of active nodes in the list L
   k                  = length(L);
   % Generate the nodes still needed for the active list of nodes
   pred               = MILP.pred(1:NODE,1);
   p                  = 1;   % Save root node as first node in the tree                     
   % Generate the search tree for each active node
   for i=1:length(L)
       j = L(i);
       while pred(j) > 1     % Generate nodes preceding list node L(i)
          j = pred(j);
          p = [j;p];
       end
   end
   % Find the unique nodes to keep, sorted, will be nodes 1,...,M in new numbering
   u             = unique([p',L]);
   M             = length(u);
   if IterPrint > 0
      fprintf('Squeeze down to  %d nodes\n',M);
   end
   % Allocate arrays
   if SOL > 0
      XX              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      HS              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      XX(:,1:M)       = MILP.XX(:,u);
      HS(:,1:M)       = MILP.HS(:,u);
   else
      XX              = spalloc(n,M+MxAlloc,min(40000,n*MaxIter));
      XX(:,1:M)       = MILP.XX(:,u);
      HS              = [];
   end
   BB                 = spalloc(n, M+MxAlloc,min(40000,n *MaxIter));
   BB(:,1:M)          = MILP.BB(:,u);

   % Lowest value possible at node, from NLP relaxation
   f_min                = Inf*ones(M+MxAlloc,1);      
   f_min(1:M,1)         = MILP.f_min(u,1);
   % Depth in tree for each node
   Depth                = zeros(M+MxAlloc,1);      
   Depth(1:M,1)         = MILP.Depth(u,1);
   % Number of integer components
   Icomp                = zeros(M+MxAlloc,1);      
   Icomp(1:M,1)         = MILP.Icomp(u,1);

   pred                 = zeros(M+MxAlloc,1);
   Problem.xVar         = zeros(M+MxAlloc,1);
   Problem.Bound        = zeros(M+MxAlloc,1);
   Problem.xBound       = zeros(M+MxAlloc,1);

   Problem.xVar(1:M,1)  = MILP.Problem.xVar(u,1);
   Problem.Bound(1:M,1) = MILP.Problem.Bound(u,1);
   Problem.xBound(1:M,1)= MILP.Problem.xBound(u,1);
   % New node numbers Lnew for active nodes in list L
   [Lval,iu,Lnew]       = intersect(L,u);
   % Find new node number for each preceding node
   p                    = MILP.pred(u,1);
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
   v_k                  = [];
   y                    = [];

   SumIter              = TotIter;
   Iter                 = 0;

elseif WarmStart == 1 
   % Get information from last run, just keep nodes 1:NODE
   xIPMin    = MILP.xIPMin;
   fIPMin    = MILP.fIPMin;
   fIPIter   = MILP.fIPIter;
   vIPMin    = MILP.vIPMin;
   L         = MILP.L;
   NODE      = MILP.NODE;
   TotIter   = MILP.TotIter;
   Gap       = MILP.Gap; % Not needed, recomputed below
   SumMinIt  = MILP.SumMinIt;
   fDual     = MILP.fDual;
   time      = MILP.time ;
   if isempty(MILP.HS)
      SOL               = 0;
   else
      SOL               = 1;
      ProbNLP.WarmStart = 1;
   end

   M             = NODE;
   u             = 1:M;
   if IterPrint > 0
      fprintf('Warm Start with %d nodes.\n',M);
   end
   % Allocate arrays
   if SOL > 0
      XX              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      HS              = spalloc(nb,M+MxAlloc,min(40000,nb*MaxIter));
      XX(:,1:M)       = MILP.XX(:,u);
      HS(:,1:M)       = MILP.HS(:,u);
   else
      XX              = spalloc(n,M+MxAlloc,min(40000,n*MaxIter));
      XX(:,1:M)       = MILP.XX(:,u);
      HS              = [];
   end
   BB                 = spalloc(n, M+MxAlloc,min(40000,n *MaxIter));
   BB(:,1:M)          = MILP.BB(:,u);

   % Lowest value possible at node, from NLP relaxation
   f_min                = Inf*ones(M+MxAlloc,1);      
   f_min(1:M,1)         = MILP.f_min(u,1);
   % Depth in tree for each node
   Depth                = zeros(M+MxAlloc,1);      
   Depth(1:M,1)         = MILP.Depth(u,1);
   % Number of integer components
   Icomp                = zeros(M+MxAlloc,1);      
   Icomp(1:M,1)         = MILP.Icomp(u,1);

   pred                 = zeros(M+MxAlloc,1);
   pred(1:M,1)          = MILP.pred(u,1);

   Problem.xVar         = zeros(M+MxAlloc,1);
   Problem.Bound        = zeros(M+MxAlloc,1);
   Problem.xBound       = zeros(M+MxAlloc,1);

   Problem.xVar(1:M,1)  = MILP.Problem.xVar(u,1);
   Problem.Bound(1:M,1) = MILP.Problem.Bound(u,1);
   Problem.xBound(1:M,1)= MILP.Problem.xBound(u,1);

   % Update fLB, lower bound on the solution
   fLB                  = min(f_min(L));

   Inform               = 0;
   ExitFlag             = 0;
   v_k                  = [];
   y                    = [];
   SumIter              = TotIter;
   Iter                 = 0;
   WarmStart            = 1;
else
   WarmStart            = 0;
   if SaveFreq >= 0
      if exist('mipSolveSave.mat')
         delete('mipSolveSave.mat');
      end
   end
   % Initial step
   if PriLev > 2 & PriLevOpt >= 1
      fprintf('=== mipSolve:    LP relaxation. ')
      fprintf('Call %s solving Phase I and II:\n',SolverLP)
   end


   ResultLP  = tomRunMini(SolverLP,ProbLP);

   ExitFlag  = ResultLP.ExitFlag;
   Inform    = ResultLP.Inform;
   Its       = ResultLP.Iter;

   if isempty(x)
      Icomp(1)     = 0;
      Ridx         = IntVars;
   else
      x            = min(x_U,max(x_L,ResultLP.x_k));
      [IC,Iidx]    = IntComp(x,IntVars,eps_I);
      Icomp(1)     = IC;
      Ridx         = IntVars(~Iidx);
   end

   v_k       = ResultLP.v_k;
   if ExitFlag==0
      y      = v_k(n+1:n+dLin);
   else
      y      = [];
   end

   if isempty(ResultLP.SOL) 
      SOL=0;
%  elseif strcmpi(SolverDLP,'minos') | strcmpi(SolverDLP,'sqopt') | ... 
%         strcmpi(SolverDLP,'lpopt') | strcmpi(SolverDLP,'lssol') | ... 
%         strcmpi(SolverDLP,'nlssol') | strcmpi(SolverDLP,'npsol') | ... 
%         strcmpi(SolverDLP,'snopt') | strcmpi(SolverDLP,'qpopt')
   else
      SOL=1;
      ProbDLP.WarmStart=1;
   end

   % Variables on upper bound must be set as basic variables for DualSolve
   B        = ResultLP.QP.B;
   fLB      = ResultLP.f_k;
   Iter     = ResultLP.Iter;
   SumMinIt = Iter;

   if ExitFlag > 0 & ExitFlag ~= 3
      if PriLev >= 2
         fprintf('No solution found to LP relaxation. ')
         fprintf('ExitFlag %d\n', ExitFlag)
      end
      x(IntVars)          = round(x(IntVars));
      ResultMIP.x_k       = x;
      ResultMIP.f_k       = fLB;
      ResultMIP.v_k       = v_k;
      ResultMIP.QP.B      = B;
      ResultMIP.QP.y      = y;
      ResultMIP.ExitFlag  = 4;
      ResultMIP.ExitText  = ExitText(4,WarmStart,1,'');
      ResultMIP.Iter      = 0;
      ResultMIP.MinorIter = SumMinIt;
      ResultMIP           = endSolveMini(Prob,ResultMIP);
      return;
   end
   % Allocate arrays
   BB    = spalloc(n, MxAlloc,min(40000,n *MaxIter));
   if SOL > 0
      XX = spalloc(nb,MxAlloc,min(40000,nb*MaxIter));
      HS = spalloc(nb,MxAlloc,min(40000,nb*MaxIter));
   else
      XX = spalloc(n,MxAlloc,min(40000,n*MaxIter));
      HS = [];
   end

   % Save solution for 1st node
   BB(1:n,1)     = sparse(B(:));
   if SOL==0
      XX(1:n,1)  = sparse(x);
   else
      XX(1:nb,1) = sparse(ResultLP.SOL.xs);
      HS(1:nb,1) = sparse(ResultLP.SOL.hs);
   end


   % Initialization

   L    = 1;             % Node List
   NODE = 1;

   % Lowest value possible at node, from LP relaxation
   f_min          = Inf*ones(MxAlloc,1);      
   f_min(1)       = fLB;
   % Depth in tree for each node, root node is 0
   Depth          = zeros(MxAlloc,1);      
   % Number of integer components
   Icomp          = zeros(MxAlloc,1);      

   fDual          = -Inf;      % Current best dual value for feasible points

   pred           = zeros(MxAlloc,1);
   Problem.xVar   = zeros(MxAlloc,1);
   Problem.Bound  = zeros(MxAlloc,1);
   Problem.xBound = zeros(MxAlloc,1);

   fDP            = -Inf;
   SumIter        = 0;
   Iter           = 0;
   fIPIter        = 0;
end

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
       ENDTREE=1;
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
      if fIPMin == Inf % Empty feasible set
         if PriLev >= 2
            fprintf('\nTotal number of LP iterations = %d\n',SumMinIt)
            disp('----------------------------------')
            disp('No feasible solution to IP problem')
            disp('----------------------------------')
         end
         ResultMIP.x_k       = [];
         ResultMIP.f_k       = fLB;
         ResultMIP.v_k       = v_k;
         ResultMIP.QP.B      = [];
         ResultMIP.QP.y      = y;
         ResultMIP.ExitFlag  = 2;
         ResultMIP.ExitText  = ExitText(2,WarmStart,SumIter+Iter,DualGap);
         ResultMIP.Iter      = max(1,Iter);
         ResultMIP.MinorIter = SumMinIt;
         ResultMIP.DualGap   = Gap;
         ResultMIP           = endSolveMini(Prob,ResultMIP);
         return;
      else
         ResultMIP.x_k       = xIPMin;
         ResultMIP.f_k       = fIPMin;
         ResultMIP.v_k       = vIPMin;
         ResultMIP.QP.B      = BIPMin;
         ResultMIP.QP.y      = yIPMin;
         if GapOK == 0
            ResultMIP.ExitFlag=0;
            ResultMIP.ExitText = ExitText(0,WarmStart,SumIter+Iter,DualGap);
         else
            ResultMIP.ExitFlag=0;
            ResultMIP.ExitText = ExitText(6,WarmStart,SumIter+Iter,DualGap);
         end
         ResultMIP.Iter      = max(1,Iter);
         ResultMIP.MinorIter = SumMinIt;
         % ResultMIP.DualGap = 0; % Gap should be 0 when converging
         ResultMIP.DualGap   = Gap;
         if IterPrint
            PrintIter(2,Iter,[],[],fIPMin,fIPIter,fLB,L,0,0,x,[],[],...
                      Inform,ExitFlag,Its,fDual,BigNum,'');
            fprintf('--- mipSolve LP Branch & Bound converged! ')
            fprintf('Its(nodes visited) = %d. ',Iter)
            fprintf('Total LP Its = %d. ',SumMinIt)
            fprintf('Optimal f(x) =%22.16f.',fIPMin);
            fprintf('\n');
            if PriLev >= 2
               if all(xIPMin == round(xIPMin))
                  xprinti(xIPMin,'x:',7,15)
               else
                  xprint(xIPMin,'x:')
               end
               if PriLev >= 4
                  xprinti(BIPMin,'B:');
               end
            end
         end
         ResultMIP           = endSolveMini(Prob,ResultMIP);
         return;
      end
   elseif GapOK == 1
      ResultMIP.x_k       = xIPMin;
      ResultMIP.f_k       = fIPMin;
      ResultMIP.v_k       = vIPMin;
      ResultMIP.QP.B      = BIPMin;
      ResultMIP.QP.y      = yIPMin;
      ResultMIP.ExitFlag  = 0;
      ResultMIP.ExitText  = ExitText(6,WarmStart,SumIter+Iter,DualGap);
      ResultMIP.Iter      = max(1,Iter);
      ResultMIP.MinorIter = SumMinIt;
      ResultMIP.DualGap   = Gap;
      if IterPrint
         PrintIter(3,Iter,[],[],fIPMin,fIPIter,fLB,L,0,0,x,[],[],...
                   Inform,ExitFlag,Its,fDual,BigNum,'');
         fprintf('--- mipSolve: User defined duality gap reached. ');
         fprintf('Its(nodes visited) = %d. ',Iter)
         fprintf('Total LP Its =%d. ',SumMinIt)
         fprintf('Optimal f(x) =%22.16f.',fIPMin);
         fprintf('\n');
         if PriLev >= 2
            if all(xIPMin == round(xIPMin))
               xprinti(xIPMin,'x:',7,15)
            else
               xprint(xIPMin,'x:')
            end
            if PriLev >= 4
               xprinti(BIPMin,'B:');
            end
         end
      end
      ResultMIP           = endSolveMini(Prob,ResultMIP);
      return;
   end
   if cputime-TIME0 > MaxCPU, cpumax = 1; break; end
   Iter = Iter+1;
   if PriLev > 2
      fprintf('\n========== mipSolve ==========   Iteration %6.0f\n',Iter);
   end
   % Problem selection and relaxation
   [i,L] = BBnode(NodeSel,L,pred,Depth,Icomp,fIPMin);

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

   if i == 1     % 1st step
      x        = full(XX(1:n,1));
      B        = full(BB( : ,1));

      Level    = 0;
      fNode    = fLB;
      f_min(1) = fNode;
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
      B = full(BB( : ,j));

      x(k) = max(x_L(k),min(x_U(k),x(k))); % Adjust optimal x for parent node inside x_L,x_U

      if PriLev > 2
         xprint(x,'x:')
      end

      if 1 & (SOL > 0)
         ProbDLP.SOL.xs    = full(XX(1:nb,j));
         ProbDLP.SOL.hs    = full(HS(1:nb,j));
      else
         ProbDLP.WarmStart = 0;
      end

      ProbDLP.P            = i; 
      ProbDLP.QP.B         = B;
      ProbDLP.x_0          = x;
      ProbDLP.x_L          = x_L; 
      ProbDLP.x_U          = x_U; 
      ProbDLP.QP.DualLimit = fIPMin;     % Try stop dual iterations early

      Result               = tomRunMini(SolverDLP,ProbDLP);

      ExitFlag             = Result.ExitFlag;
      Inform               = Result.Inform;
      if (ExitFlag == 6 | ExitFlag == 1) & strcmpi('DualSolve',SolverDLP)
         % MipSolve gave bad starting point (or too many iterations.
         % Try a standard LP solver
         Result    = tomRunMini('qld',ProbDLP);
         Inform    = Result.Inform;
      end

      if isempty(x)
         Icomp(i)     = 0;
         Ridx         = IntVars;
      else
         x            = min(x_U,max(x_L,Result.x_k));
         [IC,Iidx]    = IntComp(x,IntVars,eps_I);
         Icomp(i)     = IC;
         Ridx         = IntVars(~Iidx);
      end

      SumMinIt     = SumMinIt+Result.Iter;
      v_k          = Result.v_k;
      %y           = v_k(n+1:length(v_k));
      y            = Result.y_k;
      fDP          = Result.f_k;
      B            = Result.QP.B;
      Its          = Result.Iter;

      if ExitFlag == 0 | ExitFlag == 3
         BB(:,i) = sparse(B);
         if SOL == 0 | isempty(Result.SOL)
            XX(1:n,i) = sparse(x);
         else
            XX(:,i)   = sparse(Result.SOL.xs);
            HS(:,i)   = sparse(Result.SOL.hs);
         end

         fNode    = Result.f_k;        % fNode    = c'*x;
         f_min(i) = fNode;
         fLB      = min(f_min([i,L])); % Update best possible primal

      else
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
         else
            [IntVar,xC,xI] = BranchVar(VarSel,x,Ridx,ProbDLP.QP.c,eps_I,PriLev);
         end
      end
      if IntVar
         if fDP > fDual, fDual=fDP; end
         if fNode <  fIPMin  % Check if new integer value lower than best found
            fIPMin           = fNode;
            fIPIter          = Iter;
            x(IntVars)       = round(x(IntVars));
            xIPMin           = min(x_U,max(x_L,x(1:n)));
            BIPMin           = B(1:n);
            yIPMin           = y;
            vIPMin           = v_k;
            DelL             = f_min(L) >= fIPMin;
            XX(:,L(DelL))    = 0;
            BB(:,L(DelL))    = 0;
	    if SOL > 0
               HS(:,L(DelL))    = 0;
            end

            L = L(~DelL);

            if IterPrint
               PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),nI-Icomp(i),x,xC,xI,...
                         Inform,ExitFlag,Its,fDual,BigNum,'End of tree. Found IP: ');
            end
            if abs(fIPMin - ProbDLP.QP.c(1:n)'*xIPMin) < 1E-10
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
               GapOK = 1;
            end
         end % if fNode <  fIPMin
         ENDTREE=1;
      else % end if IntVar
         % Apply heuristics
         RndHeu = 0;
         if ROUNDH                 % Check rounding heuristic
            xRS                    = x;
            xRS(IntVars)           = round(x(IntVars));
            [LV LVi]               = consViol(ProbDLP.A*xRS,b_L,b_U,absviol);
            if LV < bTol
               xRS                 = min(x_U,max(x_L,xRS));
               fRS                 = xRS'*ProbDLP.QP.c;
               if fRS < fIPMin
                  RndHeu           = 1;
                  fIPMin           = fRS;
                  fIPIter          = Iter;
                  xIPMin           = xRS;
                  BIPMin           = B(1:n);
                  yIPMin           = [];
                  vIPMin           = [];
                  DelL             = f_min(L) >= fIPMin;
                  XX(:,L(DelL))    = 0;
                  BB(:,L(DelL))    = 0;
                  if SOL > 0
                     HS(:,L(DelL)) = 0;
                  end
                  L                = L(~DelL);
                  if IterPrint
                     PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),nI-Icomp(i),x,xC,xI,...
                               Inform,ExitFlag,Its,fDual,BigNum,'Rounding Heuristic: ');
                  end
               end
            end
         end
         if KNAPSACK > 0 & RndHeu == 0
            % Apply rounding knapsack heuristic
            xKS          = x;
            xKS(IntVars) = floor(x(IntVars)+eps_I*max(1,abs(x(IntVars))));
            if nMK > 0
               xKS(xMK) = 0;                                    % Set slacks to 0
               xKS(xMK) = b_U(bMK) - ProbDLP.A(bMK,:)*xKS;      % Find slack values
               xKS(xIK) = round(xKS(xIK));                      % Must round integer slacks
               xKS(xMK) = max(x_L(xMK),min(x_U(xMK),xKS(xMK))); % Inside bounds
            end
            [LV LVi]               = consViol(ProbDLP.A*xKS,b_L,b_U,absviol);
            if LV < bTol
               xKS                 = min(x_U,max(x_L,xKS));
               fKS                 = xKS'*ProbDLP.QP.c;
               if fKS < fIPMin
                  fIPMin           = fKS;
                  fIPIter          = Iter;
                  xIPMin           = xKS;
                  BIPMin           = B(1:n);
                  yIPMin           = [];
                  vIPMin           = [];
                  DelL             = f_min(L) >= fIPMin;
                  XX(:,L(DelL))    = 0;
                  BB(:,L(DelL))    = 0;
                  if SOL > 0
                     HS(:,L(DelL)) = 0;
                  end
                  L                = L(~DelL);
                  if IterPrint
                     PrintIter(4,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,sum(DelL),nI-Icomp(i),x,xC,xI,...
                               Inform,ExitFlag,Its,fDual,BigNum,'Knapsack Heuristic: ');
                  end
               end
            end
         end
      end
   end
   % Compute duality gap
   if ~isinf(fIPMin)
      if fIPMin ~= 0
         Gap =  max(0,fIPMin-fLB)/abs(fIPMin);
      else
         Gap =  max(0,fIPMin-fLB);
      end
      % Check if duality gap <= DualGap
      if Gap <= DualGap
         GapOK = 1;
      else
         GapOK = 0;
      end
   end
   if fNode > fIPCut             % Cut the tree, user defined
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
                Inform,ExitFlag,Its,fDual,BigNum,'');
   end
   if SaveFreq >0 & mod(Iter,SaveFreq) == 0
      time         = fix(clock);
      try
         TotIter   = SumIter+Iter;
         save('mipSolveSave.mat','Name', 'xIPMin','fIPMin','fIPIter','vIPMin',...
                            'XX','HS','BB', ...
                            'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                            'Depth','Icomp',...
                            'Gap','SumMinIt','fDual','time');
      catch
         warning('Failed to save warm start information to mipSolveSave.mat');
         err = lasterror;
         disp(err.message);
         disp('Warm start information is available in Result.mipSolve');
      end
   end
end % while
if PriLev >= 1
   fprintf('\n--- TOO MANY Branch and Bound ITERATIONS. ITER = %d\n',Iter)
   fprintf('\n--- Total number of LP iterations = %d\n',SumMinIt)
end
ResultMIP.x_k                = xIPMin;
ResultMIP.f_k                = fIPMin;
ResultMIP.v_k                = vIPMin;
ResultMIP.QP.B               = BIPMin;
ResultMIP.QP.y               = yIPMin;
% Warm start info saved
ResultMIP.mipSolve.Name      = Name;
ResultMIP.mipSolve.xIPMin    = xIPMin;
ResultMIP.mipSolve.fIPMin    = fIPMin;
ResultMIP.mipSolve.fIPIter   = fIPIter;
ResultMIP.mipSolve.vIPMin    = vIPMin;
ResultMIP.mipSolve.XX        = XX(:,1:NODE);
if SOL > 0
   ResultMIP.mipSolve.HS     = HS(:,1:NODE);
else
   ResultMIP.mipSolve.HS     = HS;
end
ResultMIP.mipSolve.BB        = BB(:,1:NODE);
ResultMIP.mipSolve.L         = L;
ResultMIP.mipSolve.NODE      = NODE;

TotIter                      = SumIter+Iter;
ResultMIP.mipSolve.TotIter   = TotIter;
ResultMIP.mipSolve.f_min     = f_min(1:NODE);
ResultMIP.mipSolve.pred      = pred(1:NODE);
ResultMIP.mipSolve.Problem.xVar     = Problem.xVar(1:NODE);
ResultMIP.mipSolve.Problem.Bound    = Problem.Bound(1:NODE);
ResultMIP.mipSolve.Problem.xBound   = Problem.xBound(1:NODE);
% ResultMIP.mipSolve.Problem = Problem;
ResultMIP.mipSolve.Depth     = Depth;
ResultMIP.mipSolve.Icomp     = Icomp;
ResultMIP.mipSolve.Gap       = Gap;
ResultMIP.mipSolve.SumMinIt  = SumMinIt;
ResultMIP.mipSolve.fDual     = fDual;
time                         = fix(clock);
ResultMIP.mipSolve.time      = time;

try
   save('mipSolveSave.mat','Name', 'xIPMin','fIPMin','fIPIter','vIPMin',...
                      'XX','HS','BB', ...
                      'L','NODE','TotIter', 'f_min','pred', 'Problem',...
                      'Depth','Icomp',...
                      'Gap','SumMinIt','fDual','time');
catch
   warning('Failed to save warm start information to mipSolveSave.mat');
   err = lasterror;
   disp(err.message);
   disp('Warm start information is available in Result.mipSolve');
end

if cpumax
   ResultMIP.ExitFlag=99;
   ResultMIP.ExitText=ExitText(9,WarmStart,TotIter,DualGap);
else
   ResultMIP.ExitFlag=1;
   ResultMIP.ExitText=ExitText(1,WarmStart,TotIter,DualGap);
end
ResultMIP.Iter      = Iter;
ResultMIP.MinorIter = SumMinIt;
ResultMIP.DualGap   = Gap;
ResultMIP           = endSolveMini(Prob,ResultMIP);

% ------------------------------
function Text = ExitText(Inform,WarmStart,Its,DualGap)
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
        Text = ['No solution found to LP relaxation. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'No solution found to LP relaxation';
     end
   case 5
     if WarmStart
        Text = ['Illegal x_0 found in LP relaxation. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = 'Illegal x_0 found in LP relaxation';
     end
   case 6
     if WarmStart
        Text = ['User defined duality gap ' num2str(DualGap) ' reached. ', ...
          'Tried in total ' num2str(Its) ' iter.'];
     else
        Text = ['User defined duality gap ' num2str(DualGap) ' reached'];
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
function [h_L1, h] = consViol(z,L,U,absviol)
% -----------------------------------------------
if absviol == 1
   h =-min(0,z-L)-min(0,U-z);
else
   h =-min(0,(z-L)./max(1,abs(L)))-min(0,(U-z)./max(1,abs(U)));
end
h_L1 = sum(h);

% -----------------------------------------------
function PrintIter(Type,Iter,Level,fNode,fIPMin,fIPIter,fLB,L,DelL,nInts,x,xC,xI,...
                   Inform,ExitFlag,Its,fDual,BigNum,Txt1);
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
         fprintf('fLP %9.0f ',fNode);
      else
         fprintf('fLP %9.3f ',fNode);
      end
   else
      fprintf('%s',blanks(14));
   end
else
   fprintf('%s%s',blanks(27-length(Txt1)),Txt1);
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
      fprintf('%8.3f%% ',100*max(0,fIPMin-fLB)/abs(fIPMin));
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
   fprintf('Its%4d',Its);

   if ~isinf(fDual) 
      if BigNum
         fprintf(' fDual %8.0f',fDual);
      else
         fprintf(' fDual %8.2f',fDual);
      end
   else
      fprintf('%s',blanks(15));
   end
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
% 981111  hkh  Change to use lpSolve (now called lpSimplex)
% 981114  hkh  Change to Prob/Result format for lpSolve
% 981115  hkh  Setting lower bounds as one, not -Inf.
% 981119  hkh  Change field to Prob.QP.B. Errors in use of B
% 981123  hkh  Change name to mipSolve. 
% 981123  hkh  Change optPar(28) to ExitFlag on 3 places. max_iter to MaxIter
% 990419  hkh  Test on idx = find(B(1:Ivar)) empty (= int vars 0, non basic).
% 990804  hkh  Change to Prob/Result input-output format
% 990810  hkh  Revised for v2.0. Using new DualSolve and lpSolve.
% 990901  hkh  Calling general SolveDLP.
% 990913  hkh  Safeguard against nan and inf in x_0
% 991222  hkh  Chcck if c is empty, stop.
% 000916  hkh  Add text about convergence in ExitText
% 010407  hkh  IntVars and Ivars made correct
% 020204  hkh  Pick up Prob.SOL.optPar/PrintFile and use if SOL solvers
% 020304  hkh  Always use MINOS instead of LPOPT for LP sub problems.
% 020304  hkh  Field error setting Prob.SOL.optPar. Print Inform from solver
% 021223  hkh  Handle lower and upper bounds on linear constraints
% 030107  hkh  Avoid allocating too much memory for XX, HS, BB
% 030107  hkh  Prob.b_L empty now gives Prob.b_U (as comments say)
% 030107  hkh  Correct comments
% 030203  hkh  Correct = b_U to <= b_U
% 030309  hkh  Ensure all x inside bounds: x = min(x_U,max(x_L,Result.x_k));
% 030823  hkh  Fix bugs in handling duality gap, incl. printout last iteration
% 030831  hkh  Duality gap not checked at start. Check feasibility of xIPMin
% 031111  hkh  Add Duality gap as output field Result.DualGap
% 031128  hkh  Gap undefined if fIPMin given. Set DualGap=0 if convergence
% 040111  hkh  Add call to inisolve and endSolve
% 040425  hkh  Add input SmallA, if 1 detect and remove small A elements
% 040425  hkh  New option: Test for max CPU Time used (cputime > Prob.MaxCPU)
% 041123  hkh  Check MIP fields to handle LP problems
% 041222  med  Safeguard added for x_0
% 050117  med  mlint revision
% 050603  hkh  Set x_0 to zero, crash if Init File provides empty x_0
% 051212  hkh  Generalized KNAPSACK>0 option to handle non-standard form
% 051212  hkh  VarWht was wrongly used instead of VarWeight on some places
% 060729  hkh  fLB was not dynamically updated, i.e. duality gap to big
% 060729  hkh  round fIPMin at some places to avoid close-0-blowup when dividing
% 060729  hkh  Define and print dual gap as nonnegative value
% 070222  hkh  Revise IntVars handling, use new format
% 070907  hkh  SolverLP/DLP picked from list, 1st with license, avoid GetSolver
% 080606  med  Switched to iniSolveMini
% 080607  hkh  Use endSolveMini, Use lp-minos not MINOS, tomRun not tomSolve
% 080607  med  Switched to tomRunMini
% 091004  hkh  Change input Prob.Solver.Alg to Prob.MIP.NodeSel, change default from 0 to 2
% 091004  hkh  Add Warm Start, and time stamp in iteration print
% 091005  hkh  Add depth in tree (Level) in iteration print
% 091007  hkh  Add printout of length(L) if IterPrint
% 091007  hkh  Handle fIP, without xIP, cutting nodes in tree search
% 091014  hkh  Major revision, similar to minlpSolve. Better NodeSel, use Depth (depth in tree)
%              and Icomp (# of components integer valued)
% 091020  hkh  Avoid catch ME, use lasterror
% 091022  hkh  Safeguard xIPMin = min(x_U,max(x_L,x(1:n)))
% 091022  hkh  Revise and fix errors in handling of input (xIP,fIP)
% 091022  hkh  Use WarmStart == 1: keep all nodes, WarmStart == 2 too costly when many nodes.
% 091022  hkh  Change default value for SaveFreq till -1, no save to file.
% 100217  hkh  Move ENDTREE=1 outside of Gap-computation, tree always ended when int.sol found
% 100217  hkh  Print value of DualGap in ExitText printout
% 100217  hkh  Use abs(fIPMin) in IterPrint relative Gap computation, avoid wrong sign
% 100217  hkh  Move Gap computation, compute every iter, also after heuristic or update of fLB.
% 100605  hkh  Bug in mod(). Must check SaveFreq >0 when doing mod(Iter,SaveFreq)
