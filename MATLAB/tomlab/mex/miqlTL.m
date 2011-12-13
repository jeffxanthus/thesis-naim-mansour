% TOMLAB MIQL MIQP/MILP Solver
%
% function Result = miqlTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj.
% PriLevOpt Print level in solver (0-4)
%           0 - No output
%           1 - Root QP and final performance summary.
%           2 - Additional branch and bound iterations counter.
%           3 - Additional total output from all generated subproblems.
%           4 - Additional problem formulations.
%
% -----------------------------------------------
% Fields used in Prob.MIP:
% -----------------------------------------------
% IntVars      Vector of indices of the integer variables.
% -----------------------------------------------
% Fields used in Prob.MIQL:
% -----------------------------------------------
% PrintFile    Name of MIQL Print file. Amount and type of printing determined
%              by PriLevOpt.
%
% eps          The user has to specify the desired final accuracy
%              (e.g. 1.0D-12). The parameter value should not be smaller 
%              than the underlying machine precision. 
% acc          The accuracy to identify integer values for integer variables.
%              If acc is less than machine precision, e.g., acc = 0, then acc
%              is set to the machine precision.
%
% branchrule   1    Maximal fractional branching.
%              2    Minimal fractional branching.
% nodestrategy 1    Best of all (large search trees).
%              2    Best of two (warmstarts, less memory for search tree.)
%              4    Depth first (good warmstarts, smallest memory, many QPs.) 
% maxnodes     Maximal number of nodes, e.g., 100,000. Should be at least
%              (2*n + 2*m + 6)^2 (default).
% maxwarmsts   Maximal number of successive warmstarts, e.g., 100, to avoid
%              numerical instabilities.
% improvebound 1    Calculate improved bounds if nodestrategy = 1.
%              0    Default
% dfdirection  1    Select direction for Depth-First according to value of 
%                   Lagrange function.
%              0    Default.
% choleskymode 0    Calculate Cholesky decomposition once and reuse it.
%              1    Calculate new Cholesky decomposition whenever warmstart
%                   is not active.
% cutprocess   Control the cutting process. Default = 0.
%              0    No cuts.
%              1    Use disjunctive cuts only.
%              2    Complemented mixed integer rounding (CMIR) cuts only.
%              3    Both disjunctive cuts and CMIR cuts.
% maxdiscuts   Maximal number of rounds for disjunctive cuts. Default 1, if
%              disjunctive cuts should be generated.
% maxCMIRcuts  Maximal number of cuts for CMIR cuts. Default 1, if CMIR cuts
%              should be generated.
% heurmode     Primal heuristic mode. Default = 0.
%              0    No primal heuristics. Default
%              1    Nearest integer.
%              2    Feasibility pump.
% reducedim    Reduced quadratic program dimension. Default = n.
% totalcuts    Total number of cutting planes.
% totalbbnodes Total number of branch and bound nodes.
%
% cuttingtime  Time to spend for generating cutting planes.
% bbtime       Time to spend for branch and bound process.
% cutgap       Gap reduced by cutting planes compared to original relaxed
%              solution.
% relaxopt     Relaxed optimal value
%
% generatecut  0       No cut generation.
%              Nonzero Cut generation was applied.
% posorthant   0       No transformation.
%              Nonzero Problem is transformed to the positive orthant.
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from miql.m (similar to TOMLAB).
% Inform   MIQL information parameter.
%          -4 = Branch and Bound root QP could not be solved.
%          -3 = Relaxed QP without feasible solution.
%          -2 = A feasible solution could not be computed by maxnode
%               subproblems.
%          -1 = A feasible solution does not exist.
%           0 = Optimal solution with unique minimizer found
%           1 = A feasible solution is found, but tree search is terminated.
%           2 = Index in intVars out of bounds.
%           3 = Length of v_k too small.
%           4 = Length of working arrays too short.
%           5 = Sizes are incorrectly set.
%           6 = Integer options for Branch and Bound are incorrectly set.
%           7 = Independent lagrangian multipliers could not be calculated,
%               working arrrays are too small, increase maxnodes.
%           8 = iOut or iPrint are incorrectly set.
%           9 = Lower variable bound (x_L) greater than upper bound (x_U).
%          11 = The continuous quadratic solver QL could not solve a quadratic
%               program after a maximal number of (40*(n+m)) iterations.
%          12 = Accuracy insufficient to attain convergence
%          13 = Internal inconsistency, division by zero
%          14 = Numerical instabilites prevent successful termination of 
%               continuous quadratic solver.
%         >90 = An input parameter was invalid
%               else, constraint # not consistent with other active.
%               The problem has no feasible solution
%
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (miql).
% SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------------------------------
%
% For a problem description, see miql.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Nov 1, 2000.  Last modified Aug 19, 2009.

function Result = miqlTL(Prob)

if nargin < 1, error('miqlTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 11; % QP solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='MIQL';
Result.SolverAlgorithm='MIQL MIQP/MILP code';

PriLev=Prob.PriLevOpt;

BIG=1E20;
[bl, bu, n, m] = defblbu(Prob, inf, 1);

% Initial checks on the inputs

x_0     = Prob.x_0(:);
if length(x_0) < n, x_0=zeros(n,1); end

% Safe-guard starting point
x_0   = max(bl(1:n),min(x_0,bu(1:n)));
x_L = bl(1:n);
x_U = bu(1:n);

% Make lower and upper bounds on linear constraints correct

if m > 0
   nA = size(Prob.A,2);
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   b_L = bl(n+1:n+m);
   b_U = bu(n+1:n+m);
else
   b_L = [];
   b_U = [];
end

% Set up the constraint matrix A

if m > 0
   ix1 = find( b_L == b_U & ~isinf(b_L) );
   ix2 = find( b_L ~= b_U & ~isinf(b_U) );
   ix3 = find( b_L ~= b_U & ~isinf(b_L) );

   A   = [full(Prob.A(ix1,:)); full(-Prob.A(ix2,:)); full(Prob.A(ix3,:))];
   b   = [-b_L(ix1); b_U(ix2); -b_L(ix3)];
   mEQ = length(ix1);
else
   A   = [];
   b   = [];
   mEQ = 0;
   ix1 = [];
   ix2 = [];
   ix3 = [];
end

H     = full(Prob.QP.F);
nnObj = size(H,1); % number of nonlinear variables

% Check if any linear part
c = full(Prob.QP.c(:));

% Determine type of problem
if isempty(c) | all(c==0)
   if isempty(H) | nnz(H) == 0
      Result.f_0=0;
   else
      Result.f_0=0.5*(x_0(1:nnObj)'*H*x_0(1:nnObj));
   end
else
   if isempty(H) | nnz(H) == 0
      Result.f_0=c(1:n)'*x_0(1:n);
   else
      Result.f_0=0.5*(x_0(1:nnObj)'*H*x_0(1:nnObj)) + c(1:n)'*x_0(1:n);
   end
end

%if isempty(Prob.Name)
%   sprintf(ProbName,'Problem %d',Prob.P);
%else
%   ProbName=Prob.Name;
%end
x_L(isinf(x_L)) = -BIG;
x_U(isinf(x_U)) =  BIG;
% Define default solver options.
MIQL        = DefPar(Prob,'MIQL',[]);
MIP         = DefPar(Prob,'MIP',[]);
PriLevOpt   = DefPar(Prob,'PriLevOpt',0);

% Integer variables
IntVars  = DefPar(MIP,'IntVars',[]);
% Logical vector for integers
IV       = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('miql: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars  = find(IV);

PrintFile   = DefPar(MIQL,'PrintFile','miql.txt');
miql_eps    = DefPar(MIQL,'eps',[]);
acc         = DefPar(MIQL,'acc',[]);
% integer options
optPar(1)   = DefPar(MIQL,'branchrule'  ,-999);
optPar(2)   = DefPar(MIQL,'nodestrategy',-999);
optPar(3)   = DefPar(MIQL,'maxnodes',    -999);
optPar(4)   = DefPar(MIQL,'maxwarmsts',  -999);
optPar(5)   = DefPar(MIQL,'improvebound',-999);
optPar(6)   = DefPar(MIQL,'dfdirection', -999);
optPar(7)   = DefPar(MIQL,'choleskymode',-999);
optPar(8)   = DefPar(MIQL,'cutprocess',  -999);
optPar(9)   = DefPar(MIQL,'maxdiscuts',  -999);
optPar(10)  = DefPar(MIQL,'maxCMIRcuts', -999);
optPar(11)  = DefPar(MIQL,'heurmode',    -999);
optPar(12)  = DefPar(MIQL,'reducedim',   -999);
optPar(13)  = DefPar(MIQL,'totalcuts',   -999);
optPar(14)  = DefPar(MIQL,'totalbbnodes',-999);
% real options
optPar(15)  = DefPar(MIQL,'cuttingtime', -999);
optPar(16)  = DefPar(MIQL,'bbtime',      -999);
optPar(17)  = DefPar(MIQL,'cutgap',      -999);
optPar(18)  = DefPar(MIQL,'relaxopt',    -999);
% logical options
optPar(19)  = DefPar(MIQL,'generatecut', -999);
optPar(20)  = DefPar(MIQL,'posorthant',  -999);

% Call solver
[x_k, Inform, Iter, iState, Ax, v, Obj] = miql ( ...
        H, A, b, mEQ, c, x_L, x_U, PriLevOpt, PrintFile, ...
        IntVars, optPar, miql_eps, acc );

% miql does not fix parameters on bounds, could be slightly off
x_k = min(x_U,max(x_L,x_k));

if m > 0
   Ax = Prob.A*x_k;
end

% Set the Lagrange multipliers according to TOMLAB standard
if ~isempty(v)
   cLamda = zeros(n+m,1);
   cLamda(n+1:n+m) = v(1:m);
   cLamda(ix1) = v(n+1:n+length(ix1));
   j = n+length(ix1)+length(ix2);
   cLamda(ix2) = v(n+length(ix1)+1:j);
   for i = 1:length(ix3)
       if v(j+i) ~= 0
          cLamda(ix3(i)) = v(j+i);
       end
   end
end

switch Inform
   case -4
      ExitFlag=4;  % Infeasible (root QP)
   case -3
      ExitFlag=4;  % Infeasible (relaxed QP)
   case -2
      ExitFlag=1;  % Infeasible and maxnodes exceeded
   case -1
      ExitFlag=4;  % Infeasible
   case 0
      ExitFlag=0;  % Convergence
   case 1
      ExitFlag=1;  % Too many iterations
   case 2
      ExitFlag=2;  % Input Errors (IntVars index out of bounds)
   case 3
      ExitFlag=10; % Input Errors (length of v_k)
   case 4
      ExitFlag=10; % Input errors (length of working arrays)
   case 5
      ExitFlag=10; % Input errors (sizes)
   case 6
      ExitFlag=10; % Input errors (options for B&B)
   case 7
      ExitFlag=3;  % lagrangian multipliers could not be calculated
   case 8
      ExitFlag=10; % Input errors (iOut or iPrint incorrectly set)
   case 9
      ExitFlag=10; % Input errors (bounds overlap)
   case 11
      ExitFlag=1;  % Too many iterations (in QL, iter > 40*(n+m))
   case 12
      ExitFlag=3;  % Insufficient accuracy for convergence
   case 13
      ExitFlag=3;  % Internal inconsistency, division by zero
   case 14
      ExitFlag=3;  % Numerical instabilities prevent successful termination
                   % of continuous quadratic solver
   otherwise
      ExitFlag=4;  % Infeasible     
end
if Inform == 2 && m > 0
   % Check if feasible
   if any(b_L - Prob.optParam.bTol > Ax) | any(b_U + Prob.optParam.bTol < Ax)
      ExitFlag = 4;
   end
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;

if ~isempty(c)
   if ~isempty(H)
      Result.g_k=H*x_k+c;
   else
      Result.g_k=c;
   end
elseif isempty(H)
   Result.g_k=[];
else
   Result.g_k=H*x_k;
end

Result.FuncEv    = 0;
Result.GradEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState, Result.bState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [], Prob.optParam.xTol, ...
         Prob.optParam.bTol, [], bl, bu, 1);

switch Inform
   case -4
      Text = 'Branch and Bound root QP could not be solved';
   case -3
      Text = 'Relaxed QP without feasible solution';
   case -2
      Text = 'No feasible solution found, MAXNODE subproblems exceeded';
   case -1
      Text = 'No feasible solution exists';
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'Too many iterations';
   case 2
      Text = 'Accuracy insufficient to attain convergence';
   case 3
      Text = 'Internal inconsistency, division by zero';
   case 4
      Text = 'Length of work array too short';
   case 5
      Text = 'An input parameter was invalid';
   case 6
      Text = 'Integer options for Branch and Bound incorrectly set';
   case 7
      Text = 'Independent Lagrangian multipliers could not be calculated, increase maxnodes';
   case 8
      Text = 'IOPT and IPRINT are incorrectly set';
   case 9
      Text = 'Lower variable bound exceed upper variable bound';
   case 11
      Text = 'Solver QL could not solve QP after a maximal number of 40*(N+M) iterations';
   case 12
      Text = 'The termination accuracy is insufficient for QL to satisfy the convergence problem';
   case 13
      Text = 'QL terminated due to an internal inconsistency, division by zero';
   case 14
      Text = 'Numerical instabilities prevent successful termination QL';
   otherwise
      if Inform > 100
         Text = ...
          str2mat(['Constraint ' num2str(Inform-100) ...
                    ' not consistent with other active.'] ...
                 ,'The problem has no feasible solution');
      else
         Text = ['Unknown Error ' num2str(Inform)];
      end
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nMIQL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('MIQL: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector iState for x and constraints = \n');
      xprinti(iState(1:min(length(iState),MAX_x)),'iState: ');
   end

   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (cLamda) = \n');
      xprinte(cLamda(1:min(length(cLamda),MAX_c)),'cLamda:');
   end
end
Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 091124 bjo  Written, based on qlTL.m
% 100409 bjo  Added missing cases for Inform values