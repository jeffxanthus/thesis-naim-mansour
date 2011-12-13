% L1LinSolve.m:
%
% Finds a linearly constrained L1 solution of a function of several variables
% with the use of any suitable TOMLAB LP solver
%
% function Result = L1LinSolve(Prob, PriLev)
%
% Minimization problem:
%
%        min  sum_i |Cx - y| + alpha*|Lx|, where y is in R^m
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%
% The L1 problem is solved in L1LinSolve by rewriting the problem as
% a linear programming problem.
% A set of m additional variables v, stored as x(n+1:m), and
% a set of m additional variables z, stored as x(n+m+1:n+2*m), is added
% If damping alpha is nonzero and nonempty:
% A set of n additional variables r, stored as x(n+2m+1:n+2m+n), and
% a set of n additional variables s, stored as x(n+2m+n+1:n+2m+2*n), is added
%
%        min    sum_i(v_i+z_i + alpha(r_i+s_i))
%         x
%        s/t   x_L <=   x              <= x_U
%                0 <=   y              <= Inf
%                0 <=   z              <= Inf
%                0 <=   r              <= Inf
%                0 <=   s              <= Inf
%              b_L <= A x              <= b_U
%                y <= C x  + v - z     <= y
%                0 <= L x  + r - s     <= 0
%
% ---------------------------------------------------------------------------
%
%
% INPUT PARAMETERS
%   Prob       Structure Prob. Prob must be defined.
%              Best is to use Prob = llsAssign(.....), if using the TQ format.
%              The problem should be created in the TOMLAB linear
%              least squares format (lls)
%
%   PriLev     The second input argument. Default == 2.
%              If PriLev == 0, L1LinSolve is silent, except for error messages.
%              If > 0, L1LinSolve prints summary output about problem
%              transformation
%              L1LinSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%              PriLev == 2 displays standard output in PrintResult.
%
%   Use Prob = probInit('name of file',problem_number'); if solving
%   a predefined problem in the Init File (IF) format.
%
%   Extra fields used in Prob:
%
%   SolverL1    Name of the TOMLAB solver. See getSolver and SolverList
%               Valid names are:
%               lpSimplex, qld
%               If TOMLAB /MINOS is installed: minos, lpopt, qpopt
%               If TOMLAB /SOL is installed: sqopt, snopt, npsol
%               If TOMLAB /Xpress is installed: xpress-mp
%               If TOMLAB /Cplex is installed: cplex
%               If TOMLAB /MINLP is installed: bqpd
%               If TOMLAB /XA is installed: xa
%               If TOMLAB /KNITRO is installed: knitro
%
%   LS.C       Linear matrix C m x n
%   LS.y       Data vector y m x 1
%   LS.damp    Damping parameter alpha, a scalar (matrix  1 x 1)
%   LS.L       Symmetric n x n damping matrix L in damping term alpha*|Lx|
%
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
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.r_k, Result.J_k,
%        Result.x_0, Result.xState, Result.cState,
%        Result.v_k, is transformed back to the original problem.
%        The output in Result.Prob is the result after L1LinSolve
%        transformed the problem, i.e. the altered Prob structure
%        Therefore PrintResult cannot be used

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Dec 16, 2002.   Last modified Aug 14, 2006.

function Result = L1LinSolve(Prob, PriLev)

if nargin < 2
   PriLev = [];
   if nargin < 1
      error('L1LinSolve needs input structure Prob');
   end
end
if isempty(PriLev), PriLev = 2; end

if isfield(Prob,'SolverL1')
   Solver=Prob.SolverL1;
   Solver=deblank(Solver);
else
   Solver = [];
end

solvType=checkType('lp');

if isempty(Solver)
   Solver=GetSolver('lp',Prob.LargeScale,0); 
end

probType = solvType;
Prob.probType = probType;

n = Prob.N;
[mJ1,mJ2] = size(Prob.LS.C);
if n == 0
   n = mJ2;
end
if n > 100,  Prob.LargeScale = 1; end
if n == 0
   error('Empty matrix in linear least squares problem');
end

% Check type of damping, if any

alpha = Prob.LS.damp;

if isempty(alpha) | alpha <= 0
   DAMP = 0;
   n2   = 0;
elseif isempty(Prob.LS.L)
   DAMP = 1;
   n2   = 2*n;
else
   L = Prob.LS.L;
   if size(L,1) ~= n | size(L,2) ~= n
      error('Illegal size of Prob.LS.L, should be n x n')
   end
   if ~all(all(L==L'))
      error('Prob.LS.L is not symmetric')
   end
   DAMP = 2;
   n2   = 2*n;
end
% Safe-guard starting point
Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

% Check the matrix and vector supplied
m         = mJ1;
Prob.L1.m = m;
mm        = m+m;
if m == 0
   fprintf('Number of residuals supplied is 0\n');
   error('Illegal size of initial residual vector');
end
y = Prob.LS.y;
mmY = length(y);
if mmY == 0
   y = zeros(m,1);
elseif mmY ~= m
   fprintf('Number of residual elements is %d\n',m);
   fprintf('The number of observations in Prob.LS.y is %d\n',mmY);
   error('Should be equal, illegal size of vector Prob.LS.y');
else
   y = y(:);
end

if ~isempty(Prob.A)
   mA = size(Prob.A,1);
   Prob.A = [sparse(Prob.A),sparse(mA,mm)];
else
   mA = 0;
end

if mJ1 > 0 | mJ2 > 0
   if mJ1 ~= m
      fprintf('Number of residual elements is %d\n',m);
      fprintf('There are %d rows in Prob.LS.C\n',mJ1);
      error('Should be equal, illegal size of Prob.LS.C');
   end
   if mJ2 ~= n
      fprintf('Number of variables is %d\n',n);
      fprintf('Prob.LS.C has %d columns\n',mJ2);
      error('Should be equal, illegal size of Prob.LS.C');
   end
end

% Reformulate problem
r  = Prob.LS.C*Prob.x_0 - y;
r1 = zeros(m,1);
r2 = zeros(m,1);
ix = find(r < 0);
r1(ix) = -r(ix);
ix = find(r > 0);
r2(ix) = r(ix);
mm = 2*m;

% Add 2*m extra variables (and 2*n if damping)
Prob.N    = n + mm + n2;
Prob.x_L = [Prob.x_L;zeros(mm+n2,1)];
Prob.x_U = [Prob.x_U;Inf*ones(mm+n2,1)];
% Create initial values and linear objective
Prob.QP.F = [];
if DAMP == 0
   Prob.x_0 = [Prob.x_0;r1;r2];
   Prob.QP.c = [zeros(n,1);ones(mm,1)];
else
   Prob.x_0 = [Prob.x_0;r1;r2; max(0,Prob.x_0);min(0,Prob.x_0)];
   Prob.QP.c = [zeros(n,1);ones(mm,1);alpha(1)*ones(2*n,1)];
end
% Expand with linear constraints
if DAMP == 0
   Prob.A    = [sparse(Prob.A);sparse(Prob.LS.C),-speye(m,m),speye(m,m)];
   Prob.b_L  = [Prob.b_L;y];
   Prob.b_U  = [Prob.b_U;y];
elseif DAMP == 1
   Prob.A    = [sparse(Prob.A),sparse(mA,n2);
      sparse(Prob.LS.C),-speye(m,m),speye(m,m),sparse(m,n2); ...
      speye(n,n),sparse(n,mm),-speye(n,n),speye(n,n)];
   Prob.b_L  = [Prob.b_L;y;zeros(n,1)];
   Prob.b_U  = [Prob.b_U;y;zeros(n,1)];
else
   Prob.A    = [sparse(Prob.A),sparse(mA,n2);...
      sparse(Prob.LS.C),-speye(m,m),speye(m,m),sparse(m,n2); ...
      sparse(L),sparse(n,mm),-speye(n,n),speye(n,n)];
   Prob.b_L  = [Prob.b_L;y;zeros(n,1)];
   Prob.b_U  = [Prob.b_U;y;zeros(n,1)];
end
Prob.mLin = size(Prob.A,1);

if PriLev > 0
   fprintf('The problem has %d variables (Columns),',n);
   fprintf(' L1LinSolve adds non-negative variables (Columns) %d to %d\n',...
      n+1,n+mm);
   fprintf('The sum of the extra variables are residual functions\n');
   fprintf('\n');
   if DAMP > 0
      fprintf(' L1LinSolve adds non-negative variables (Columns) %d to %d\n',...
         n+mm+1,n+mm+n2);
      fprintf('The sum of the extra variables are the sum of x\n');
      fprintf('\n');
   end
   if mA == 0
      fprintf('The problem has %d linear constraints (Rows)\n',mA);
   else
      fprintf('The problem has %d linear constraints',mA);
      fprintf(' defined as constraint (Row) %d to %d\n',1,mA);
   end
   fprintf('\n');
end

% Remake Prob.FUNCS into an LP problem
Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');
Prob.CHECK = 0;
Result=tomRun(Solver,Prob,PriLev);

Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation
% Return only original number of variables, remove f(x) value slack
if ~isempty(Result.x_0)
   Result.x_0     = Result.x_0(1:n);
end
if ~isempty(Result.xState)
   Result.xState  = Result.xState(1:n);
end

if ~isempty(Result.v_k)
   Result.v_k  = Result.v_k(1:n);
end
if ~isempty(Result.x_k)
   Result.r_k  = Result.x_k(n+m+1:n+mm)-Result.x_k(n+1:n+m);
   Result.x_k  = Result.x_k(1:n);
end

% MODIFICATION LOG:
%
% 030524  hkh  Use tomFiles to set Prob.FUNCS
% 040125  hkh  Set Prob.mLin after changing Prob.A
% 040728  med  tomFiles used instead
% 041123  hkh  Change call to tomRun
% 041201  hkh  Removed Prob = iniSolve(Prob,8,0,0);, not a valid call
% 041201  hkh  Removed unnecessary call to ProbCheck, done by tomRun
% 041201  hkh  Force a check in ProbCheck, setting Prob.CHECK=0
% 041202  med  Added isempty check for x_0
% 041202  hkh  Check size and symmetry of Prob.LS.L
% 041202  hkh  Change comments for Prob.LS.L and Prob.LS.damp
% 041202  hkh  Define and use LS=Prob.LS, set Prob.LS = [];
% 041222  med  Safeguard added for x_0
% 050412  hkh  Avoid intermediate LS structure
% 050421  hkh  Comments on some unnecessary tests
% 050422  hkh  Avoid Lagrange multipliers in PrintResult
% 050831  med  Removed size checks on A, b_L, b_U, x_L, x_U, x_0
% 050831  med  Removed all variables not used in the code (mlint)
% 050831  med  All eye statements replaced by speye
% 050901  med  Making sure that Prob.A is sparse at all times
% 050901  med  All inputs to Prob.A are now sparse
% 060814  med  FUNCS used for callbacks instead