% linRatSolve.m:
%
% Finds a linearly constrained solution of a function of the ratio
% of two linear functions. Binary and integer variables are not supported.
%
% function Result = linRatSolve(Prob, PriLev)
%
% Minimization problem:
%
%        min  (c1 * x) / (c2 * x), where c1/c2 are in R^n
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%
% The linear ratio problem is solved by linRatSolve by rewriting the
% problem as a linear constrained optimization problem.
% n+1 variables z1 and z2(2:n+1) are needed, stored as x(1:n+1). The n original
% variables are removed so one more variable exists in the final problem.
%
%              z1 = 1/sum(c2 * x)
%              z2 = x.*z1
%
% The problem then becomes:
%
%              z1*(c1 * x) = (c1 * z1 * x) = c1 * z2
%
%        min    c1 * z2,    where c1 * z2 is in R^n (or R^n+1)
%         x
%        s/t  z1_L <=   z1             <= Inf
%                1 <=   c2.*z2         <= 1
%                0 <=   A z2 - z1*beq  <= 0  % Equalities
%             -Inf <=   A z2 - z1*b_U  <= 0  % Inequality from upper bounds
%             -Inf <=  -A z2 + z1*b_L  <= 0  % Inequality from lower bounds
%
%                0 <=   A1 z2 - z1*xeq <= 0  % Equality from bounds
%             -Inf <=   A1 z2 - z1*x_U <= 0  % Inequality from upper bounds
%             -Inf <=  -A1 z2 + z1*x_L <= 0  % Inequality from lower bounds
%
% OBSERVE the denominator c2 * x must always be positive. It is normally a
% good a idea to run the problem with both signs (multiply each side by -1).
%
% ---------------------------------------------------------------------------
%
%
% INPUT PARAMETERS
%   Prob       Structure Prob. Prob must be defined.
%              Best is to use Prob = lpAssign(.....), if using the TQ format.
%              Prob.QP.c1/c2 matrices should then be set (Prob.QP.c ignored).
%
%   PriLev     The second input argument. Default == 2.
%              If PriLev == 0, linRatSolve is silent, except for error messages.
%              If > 0, linRatSolve prints summary output about problem
%              transformation.
%              linRatSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%
%   Extra fields used in Prob:
%
%   SolverRat   Name of the TOMLAB solver. Valid names are:
%               cplex, minos, snopt, xa and more. See SolverList('lp');
%
%   QP.c1      The numerator in the objective.
%   QP.c2      The denominator in the objective.
%
%   z1_L       Lower bound for z1 (default 1e-5).
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
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.f_k
%        Result.x_0, Result.xState, Result.bState, Result.v_k, is
%        transformed back to the original problem.
%        The output in Result.Prob is the result after linRatSolve
%        transformed the problem, i.e. the altered Prob structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 5, 2005.   Last modified Oct 6, 2005.

function Result = linRatSolve(Prob, PriLev)

if nargin < 2
   PriLev = [];
   if nargin < 1
      error('linRatSolve needs input structure Prob');
   end
end

if isempty(PriLev), PriLev = 2; end

if isfield(Prob.QP,'c1')
   if isempty(Prob.QP.c1)
      error('No linear objectives given in Prob.QP.c1');
   else
      c1 = full(Prob.QP.c1(:));
   end
end

if isfield(Prob.QP,'c2')
   if isempty(Prob.QP.c2)
      error('No linear objectives given in Prob.QP.c2');
   else
      c2 = full(Prob.QP.c2(:));
   end
end

if isfield(Prob.MIP,'IntVars')
    Prob.MIP.IntVars = [];
end

if length(c2) ~= Prob.N | length(c1) ~= Prob.N
   error('Prob.QP.c1 or Prob.QP.c2 do not have the correct length');
end

if isfield(Prob,'SolverRat')
   Solver=deblank(Prob.SolverRat);
else
   Solver = [];
end

solvType=checkType('lp');
if isempty(Solver), Solver=GetSolver('lp',Prob.LargeScale,0); end

Prob.CHECK = 0; % Force a test in ProbCheck
Prob=ProbCheck(Prob,Solver,solvType);

n = Prob.N;

if isempty(Prob.x_0)
   Prob.x_0 = zeros(n,1);
end
% Safe-guard starting point
Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

% Modifying A, first check which constraints are equalities
mA = size(Prob.A,1);
% Find and remove equalities
idxEq  = find(abs(Prob.b_L-Prob.b_U) < 1e-8);
beq    = Prob.b_L(idxEq);
Prob.b_L(idxEq) = [];
Prob.b_U(idxEq) = [];
Aeq    = Prob.A(idxEq,:);
Prob.A(idxEq,:)   = [];

% Find all inequalities
idxLow = find(Prob.b_L < 1e20 & Prob.b_L > -1e20);
idxUpp = find(Prob.b_U < 1e20 & Prob.b_U > -1e20);
Alow   = Prob.A(idxLow,:);
Aupp   = Prob.A(idxUpp,:);
blow   = Prob.b_L(idxLow);
bupp   = Prob.b_U(idxUpp);

Prob.A   = [sparse([0,c2']);sparse([-beq(:),Aeq]);sparse([-bupp(:),Aupp]);sparse([blow(:),-Alow])];
Prob.b_L = [1;zeros(length(beq),1);-inf*ones(length(blow)+length(bupp),1)];
Prob.b_U = [1;zeros(length(beq),1);zeros(length(blow)+length(bupp),1)];
Prob.N = n+1;
Prob.QP.c = [0;c1];

% Adding lower/upper bounds into linear constraints
% x_L <= speye(Prob.N) <= x_U
A1     = speye(n);

% Find all double sided and single-sided zeros.
idxFix = find(abs(Prob.x_L) < 1e-8 & abs(Prob.x_U) < 1e-8)';
idxx_L = find(abs(Prob.x_L) < 1e-8 & ~(abs(Prob.x_U) < 1e-8))';
idxx_U = find(~(abs(Prob.x_L) < 1e-8) & abs(Prob.x_U) < 1e-8)';

% Find and remove equalities
idxEq  = find(abs(Prob.x_L-Prob.x_U) < 1e-8 & ~(abs(Prob.x_L) < 1e-8 & abs(Prob.x_U) < 1e-8));
xeq    = Prob.x_L(idxEq);
Prob.x_L(idxEq) = [];
Prob.x_U(idxEq) = [];
A1eq   = A1(idxEq,:);
A1(idxEq,:)   = [];

% Find all inequalities
idxLow = find(Prob.x_L < 1e20 & Prob.x_L > -1e20 & ~(abs(Prob.x_L) < 1e-8));
idxUpp = find(Prob.x_U < 1e20 & Prob.x_U > -1e20 & ~(abs(Prob.x_U) < 1e-8));
A1low   = A1(idxLow,:);
A1upp   = A1(idxUpp,:);
xlow   = Prob.x_L(idxLow);
xupp   = Prob.x_U(idxUpp);

Prob.A   = [Prob.A;sparse([-xeq(:),A1eq]);sparse([-xupp(:),A1upp]);sparse([xlow(:),-A1low])];
Prob.b_L = [Prob.b_L;zeros(length(xeq),1);-inf*ones(length(xlow)+length(xupp),1)];
Prob.b_U = [Prob.b_U;zeros(length(xeq),1);zeros(length(xlow)+length(xupp),1)];
Prob.mLin = size(Prob.A,1);

Prob.x_0 = ones(n+1,1);
Prob.x_L = -inf*ones(n+1,1);
if isfield(Prob,'z1_L')
   Prob.x_L(1,1) = Prob.z1_L;
else
   Prob.x_L(1,1) = 1e-5;
end
Prob.x_U = inf*ones(n+1,1);

Prob.x_L(idxFix+1) = 0;
Prob.x_U(idxFix+1) = 0;
Prob.x_L(idxx_L+1) = 0;
Prob.x_U(idxx_U+1) = 0;

if PriLev > 0
   fprintf('The problem has %d variables (Columns),',n);
   fprintf(' linRatSolve replaces these variable with z1 and z2(1:n)\n');
   if mA == 0
      fprintf('The problem has %d linear constraints (Rows)\n',mA);
   else
      fprintf('The problem has %d linear constraints',mA);
      fprintf(' defined as constraint (Row) %d to %d\n',1,mA);
   end
end

Prob = tomFiles(Prob,'lp_f','lp_g','lp_H');

Prob.NumDiff          = 0;
Prob.ConsDiff         = 0;
Result=tomRun(Solver,Prob,PriLev);

Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation

% Return only original number of variables, remove f(x) value slack

Result.xState  = Result.xState(1:n);
Result.bState  = Result.bState(1:mA);

if ~isempty(Result.v_k)
   Result.v_k  = Result.v_k(1:n);
end

if ~isempty(Result.x_k)
   Result.x_k  = Result.x_k(2:n+1)/Result.x_k(1);
end

Result.H_k  = [];

% MODIFICATION LOG:
%
% 051005  med  Written.
% 051006  med  Transfer of zeros in bounds to new problem
% 051006  med  General help/code review
% 080514  med  Prob.MIP.IntVars automatically purged