% infLinSolve.m:
%
% Finds a linearly constrained minimax solution of a function of several
% variables with the use of any suitable TOMLAB solver. The decision
% variables may be binary or integer.
%
% function Result = infLinSolve(Prob, PriLev)
%
% Minimization problem:
%
%        min  max D*x, where D is in R^m*n
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%
% The different objectives are stored in D row-wise.
%
% x(IntVars) are integer values, IntVars is an index set, a subset of [1:n].
%
% The linear minimax problem is solved in infLinSolve by rewriting the
% problem as a linear constrained optimization problem.
% One additional variable z, stored as x(n+1), is added
%
%        min    z
%         x
%        s/t   x_L <=   x(1:n)         <= x_U
%             -Inf <=   z              <= Inf
%              b_L <= A x              <= b_U
%             -Inf <= D*x - z*e        <= 0,  e is in R^m, e(i)=1 for all i
%
% To handle cases where a row in D*x is taken the absolute
% value of:  min max |D*x|, expand the problem with extra residuals
% with the opposite sign: [D*x; -D*x]
%
% ---------------------------------------------------------------------------
%
%
% INPUT PARAMETERS
%   Prob       Structure Prob. Prob must be defined.
%              Best is to use Prob = lp/mipAssign(.....), if using the TQ format.
%              Prob.QP.D matrix should then be set to the rows (Prob.QP.c
%              ignored).
%
%   PriLev     The second input argument. Default == 2.
%              If PriLev == 0, infLinSolve is silent, except for error messages.
%              If > 0, infLinSolve prints summary output about problem
%              transformation.
%              infLinSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%
%   Extra fields used in Prob:
%
%   SolverInf   Name of the TOMLAB solver. Valid names are:
%               cplex, minos, snopt, xa and more. See SolverList('lp');
%               or SolverList('mip');
%
%   QP.D        The rows with the different objectives.
%
%   f_Low       A lower bound on the optimal function value.
%               Not crucial, if not set default -1E300 is used
%   f_Upp       An upper bound on the optimal function value.
%               Not crucial, if not set default 1E300 is used
%
%   The rest of the fields in Prob should be defined as wanted by the
%   selected solver. See the help for the solver.
%
%   In particular:
%
%   x_0:     Starting point
%   x_L:     Lower bounds for x
%   x_U:     Upper bounds for x
%   b_L:     Lower bounds for linear constraints
%   b_U:     Upper bounds for linear constraints
%   A:       Linear constraint matrix
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%           Defines integer optimization parameters. Fields used:
%  IntVars:  
%           If empty, all variables are assumed non-integer 
%           If islogical(IntVars) (=all elements are 0/1), then
%           1 = integer variable, 0 = continuous variable.
%           If any element >1, IntVars is the indices for integer variables
%
%
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.f_k
%        Result.x_0, Result.xState, Result.bState, Result.v_k, is
%        transformed back to the original problem.
%        The output in Result.Prob is the result after infLinSolve
%        transformed the problem, i.e. the altered Prob structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Oct 4, 2005.   Last modified Feb 22, 2007.

function Result = infLinSolve(Prob, PriLev)

if nargin < 2
   PriLev = [];
   if nargin < 1
      error('infLinSolve needs input structure Prob');
   end
end

if isempty(PriLev), PriLev = 2; end

if isfield(Prob.QP,'D')
   if isempty(Prob.QP.D)
      error('No linear objectives given in Prob.QP.D');
   end
end

solvType = [];

if isfield(Prob,'SolverInf')
   Solver=deblank(Prob.SolverInf);
else
   Solver = [];
   % Use MILP solver if integer variables
   IntVars  = DefPar(Prob.MIP,'IntVars',[]);

   % Logical vector for integers
   IV = zeros(n,1);

   if isempty(IntVars)
      % No binary variables B or integer variables of type I
   elseif any(IntVars==0) | all(IntVars==1)
      % Assume binary logical vector given
      IV(1:length(IntVars)) = logical(IntVars);
   else
      if any(IntVars < 1 | IntVars > n)
         error('infLinSolve: Illegal IntVars vector');
      end
      IV(IntVars)=1;
   end
   IntVars = find(IV);
   Prob.MIP.IntVars = IntVars; % Now works when 1 more variable used
   if ~isempty(IntVars)
      Solver   = GetSolver('mip',Prob.LargeScale,0);
      solvType = checkType('mip');
   end
end

if isempty(Solver), Solver=GetSolver('lp',Prob.LargeScale,0); end
if isempty(solvType), solvType=checkType('lp'); end

Prob.CHECK = 0; % Force a test in ProbCheck
Prob=ProbCheck(Prob,Solver,solvType);

n = Prob.N;

if isempty(Prob.x_0)
   Prob.x_0 = zeros(n,1);
end

% Safe-guard starting point
Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

% Check the Prob.QP.D matrix supplied
[md,nd] = size(Prob.QP.D);

if nd ~= n
   error('Illegal number of columns in Prob.QP.D');
end

% Modifying A to include one more column
if ~isempty(Prob.A)
   mA = size(Prob.A,1);
   Prob.A = [sparse(Prob.A),sparse(mA,1)];
else
   mA = 0;
end

% Defining the new objective
Prob.QP.c = zeros(n+1,1);
Prob.QP.c(end,1) = 1;

% Modifying A to include the new constraints
Prob.A    = [sparse(Prob.A);sparse(Prob.QP.D),sparse(-1*ones(md,1))];
Prob.mLin = Prob.mLin + md;
Prob.b_L  = [Prob.b_L;-inf*ones(md,1)];
Prob.b_U  = [Prob.b_U;zeros(md,1)];
% Add one extra variables
% Largest residual based on starting point used
Prob.x_0 = [Prob.x_0;max(Prob.QP.D*Prob.x_0)];

Prob.x_L = [Prob.x_L;Prob.f_Low];

if isfield(Prob,'f_Upp')
   if isempty(Prob.f_Upp)
      Prob.x_U = [Prob.x_U;1E300];
   else
      Prob.x_U = [Prob.x_U;max(Prob.x_L(end),Prob.f_Upp)];
   end
else
   Prob.x_U = [Prob.x_U;Inf];
end
Prob.N   = n + 1;

if PriLev > 0
   fprintf('The problem has %d variables (Columns),',n);
   fprintf(' infLinSolve adds unbounded variable z as (Column) %d\n',Prob.N);
   fprintf('This extra variable z is the objective function value, where');
   fprintf(' %d <= z <= %d\n',Prob.x_L(end),Prob.x_U(end));
   fprintf('\n');
   if mA == 0
      fprintf('The problem has %d linear constraints (Rows)\n',mA);
   else
      fprintf('The problem has %d linear constraints',mA);
      fprintf(' defined as constraint (Row) %d to %d\n',1,mA);
   end
   fprintf('The problem has %d residuals, by infLinSolve ',md);
   fprintf('defined as constraint (Row) %d to %d',mA+1,mA+md);
   fprintf('\n');
end

Prob = tomFiles(Prob,'lp_f','lp_g','lp_H');

Prob.NumDiff          = 0;
Prob.ConsDiff         = 0;
Result=tomRun(Solver,Prob,PriLev);

Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation

% Return only original number of variables, remove f(x) value slack
if ~isempty(Result.x_0)
   Result.x_0     = Result.x_0(1:n);
end
Result.xState  = Result.xState(1:n);
Result.bState  = Result.bState(1:mA);

if ~isempty(Result.v_k)
   Result.v_k  = Result.v_k(1:n);
end
if ~isempty(Result.x_k)
   Result.x_k  = Result.x_k(1:n);
end

Result.H_k  = [];

% MODIFICATION LOG:
%
% 051004  med  Written.
% 051006  med  Final code review and help updates
% 060818  hkh  Use f_Low directly in Prob.x_L = [Prob.x_L;Prob.f_Low];
% 060818  hkh  Add comments about f_Low, f_Upp (set default 1E300)
% 070222  hkh  Revise IntVars handling, use new format.
