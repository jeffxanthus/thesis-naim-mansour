% TOMLAB QL QP/LP Solver
%
% function Result = qlTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj.
% PriLevOpt Print level in solver (0-1)
%           0 - No output
%           1 - Only final error message
%
% -----------------------------------------------
% Fields used in Prob.QL:
% -----------------------------------------------
% eps          The user has to specify the desired final accuracy
%              (e.g. 1.0D-12). The parameter value should not be smaller 
%              than the underlying machine precision. 
% PrintFile    Name of QL Print file. Amount and type of printing determined
%              by PriLevOpt.
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
% ExitFlag Exit status from ql.m (similar to TOMLAB).
% Inform   QL information parameter.
%           0 = Optimal solution with unique minimizer found
%           1 = Too many iterations
%           2 = Accuracy insufficient to attain convergence
%           3 = Internal inconsistency, division by zero
%           5 = An input parameter was invalid
%           else, constraint # not consistent with other active.
%           The problem has no feasible solution
%
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (ql).
% SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------------------------------
%
% For a problem description, see ql.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Nov 1, 2000.  Last modified Aug 19, 2009.

function Result = qlTL(Prob)

if nargin < 1, error('qlTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 2; % QP solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='QL';
Result.SolverAlgorithm='QL QP/LP code';

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
QL        = DefPar(Prob,'QL',[]);
PriLevOpt = DefPar(Prob,'PriLevOpt',0);
PrintFile = DefPar(QL,'PrintFile','ql.txt');
ql_eps    = DefPar(QL,'eps',[]);

% Call solver
[x_k, Inform, Iter, iState, Ax, v, Obj] = ql ( ...
        H, A, b, mEQ, c, x_L, x_U, PriLevOpt, PrintFile, ql_eps);

% ql does not fix parameters on bounds, could be slightly off
x_k = min(x_U,max(x_L,x_k));

if m > 0
   Ax = Prob.A*x_k;
end

% Get the Lagrange multipliers back to the correct place in TOMLAB
if ~isempty(v)
   cLamda = zeros(n+m,1);
   cLamda(1:n) = v(1:n);
   cLamda(n+ix1) = v(n+1:n+length(ix1));
   j = n+length(ix1)+length(ix2);
   cLamda(n+ix2) = v(n+length(ix1)+1:j);
   for i = 1:length(ix3)
       if v(j+i) ~= 0
          cLamda(n+ix3(i)) = v(j+i);
       end
   end
end

switch Inform
   case 0
     ExitFlag=0;  % Convergence
   case 1
     ExitFlag=1;  % Too many iterations
   case {2,3}
     ExitFlag=3;  % Rank problem
   case 5
     ExitFlag=10; % Input errors
%  case -1
%    ExitFlag=2;  % Unbounded
   otherwise
     ExitFlag=4;  % Infeasible
end
if Inform == 2 & m > 0
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
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'Too many iterations';
   case 2
      Text = 'Accuracy insufficient to attain convergence';
   case 3
      Test = 'Internal inconsistency, division by zero';
   case 5
      Text = 'An input parameter was invalid';
   otherwise
      Text = ...
       str2mat(['Constraint ' num2str(Inform-100) ...
                 ' not consistent with other active.'] ...
              ,'The problem has no feasible solution');
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nQL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('QL: Inform = %2d, ',Inform)
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
% 001101 hkh  Written
% 020821 hkh  Remove code appearing twice
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040102 hkh  Return only Iter ~= 0, FuncEv=GradEv=ConstrEv=0
% 040224 hkh  Check if feasible when returning Inform = 2, otherwise ExitFlag=4
% 040602 med  Help fix PriLev to PriLevOpt
% 041202 hkh  Revise calls to defblbu and StateDef, use output from defblbu
% 041222 med  Safeguard added for x_0
% 051216 med  Help updated with Inform information
% 060818 med  isnan checks removed for x_0
% 070621 frhe all(H == 0) changed to nnz(H) == 0 to avoid vector output
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
% 090605 bjo  Name changed from qld to ql
% 090605 bjo  New Inform values, 3 and > 100
% 090717 med  f_0 calculation updated
