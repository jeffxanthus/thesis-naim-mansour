% TOMLAB LSEI Least Squares Solver
%
% function Result = lseiTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% LS.C      Linear matrix m x n.
% LS.y      Data vector m x 1.
% PriLevOpt Print Level.
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% r_k      Residual vector.
% J_k      Jacobian, is just the Prob.LS.C matrix.
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial  solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status (similar to TOMLAB).
% Inform   LSEI information parameter.
% Iter     Number of iterations
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% Solver   Name of the solver (lsei).
% SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------------------------------
%
% For a problem description, see lsei.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 8, 2000.     Last modified July 22, 2011.

function Result = lseiTL(Prob)

%#function ls_f ls_g lls_H lls_r lls_J

if nargin < 1, error('lseiTL needs the Prob structure as input');end

global MAX_x % Max number of variables/constraints/resids to print
Prob.solvType = 5; % LLS solver

Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='LSEI';
Result.SolverAlgorithm='LSEI LS code';

PriLev=Prob.PriLevOpt;

% Initial checks on the inputs

[bl,bu,n,nnLin] = defblbu(Prob,Inf,1);

% Linear least squares
C = full(Prob.LS.C);
y = full(Prob.LS.y);

% The Tomlab weighting is applied in llsAssign

x_0     = Prob.x_0(:);
if length(x_0) < n, x_0=zeros(n,1); end

% Safe-guard starting point
x_0    = max(bl(1:n),min(x_0,bu(1:n)));

% Make lower and upper bounds on x correct
x_L = bl(1:n);
x_U = bu(1:n);
% Make lower and upper bounds on linear constraints correct
b_L = bl(n+1:n+nnLin);
b_U = bu(n+1:n+nnLin);

[m,nA] = size(Prob.A);

if m > 0
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
end

% Set up the constraint matrices E x = f and G x >= h

if m > 0
   ix1 = find( b_L == b_U & ~isinf(b_L) );
   ix2 = find( b_L ~= b_U & ~isinf(b_U) );
   ix3 = find( b_L ~= b_U & ~isinf(b_L) );

   E   = full(Prob.A(ix1,:));
   f   = b_L(ix1);
   G   = [full(-Prob.A(ix2,:)); full(Prob.A(ix3,:))];
   h   = [-b_U(ix2); b_L(ix3)];
else
   E   = [];
   f   = [];
   G   = [];
   h   = [];
end
ixL = find(~isinf(x_L));
ixU = find(~isinf(x_U));

mL = length(ixL);
if mL > 0
   GL = sparse(1:mL,ixL,ones(mL,1),mL,n);
   G  = [G;GL];
   h  = [h;x_L(ixL)];
end
mU = length(ixU);
if mU > 0
   GU = sparse(1:mU,ixU,-ones(mU,1),mU,n);
   G  = [G;GU];
   h  = [h;-x_U(ixU)];
end
G = full(G);

% Linear least squares
r = C*x_0-y;
Result.f_0=0.5*(r'*r);

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('qpopt',2,...
%         nnObj, 0, size(Prob.A,1), Prob);

epsRank  = Prob.optParam.eps_Rank;
% epsEQ    = Prob.optParam.bTol;
CompCov  = 0;
ScaleCov = 1;
ColScale = 1;
D        = [];
epsEQ    = [];

[x_k, rNormE, rNormL, Inform, Cov, eqRank, rlsRank, Iter] = lsei ( ...
    C, y, E, f, G, h, CompCov, ScaleCov, ColScale, D, epsEQ, epsRank);

if m > 0 & Inform == 0
   Ax = Prob.A*x_k;
else
   Ax = zeros(m,1);
end

switch Inform
   case {0,1}
     ExitFlag=0;  % Convergence
%  case 1
%    ExitFlag=1;  % Too many iterations
%  case {2}
%    ExitFlag=3;  % Rank problem
   case {4}
     ExitFlag=10; % Input errors
%  case -1
%    ExitFlag=2;  % Unbounded
   case {2,3}
     ExitFlag=4;  % Infeasible
   otherwise
     ExitFlag=999;
end
% Linear least squares
r = C*x_k-y;
Result.r_k = r;
Result.f_k = 0.5*(r'*r);
Result.J_k = C;

Result.x_k = x_k;
Result.x_0 = x_0;
%Result.v_k = cLamda;

if ~isempty(C)
   Result.g_k=C'*r;
else
   Result.g_k=zeros(n,1);
end

Result.FuncEv    = Iter;
global n_r n_J
n_r              = 2;
n_J              = 2;
Result.ResEv     = 2;
Result.JacEv     = 2;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState, Result.bState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [], Prob.optParam.xTol, ...
         Prob.optParam.bTol, [], bl, bu,1);

switch Inform
   case 0
      Text = ...
       str2mat('Both equality and inequality constraints are compatible',...
                'and have been satisfied.');
   case 1
      Text = ...
       str2mat('Equality constraints are contradictory.',...
                'A generalized inverse solution of Ex=f was used', ...
                'to minimize the residual vector length f-Ex.');
   case 2
      Text = 'Inequality constraints are contradictory.';
   case 3
      Text = 'Both equality and inequality constraints are contradictory.';
   case 4
      Text = 'Error in input parameters';
   otherwise
      Text = 'lsei: System error, illegal return parameter';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLSEI solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LSEI: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');


   fprintf('Objective function at solution x %36.18f\n\n',Result.f_k);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:  ');
   end
end

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 030508 hkh Weight problem according to Tomlab standard using WeightType
% 040102 hkh Revision for v4.2, call iniSolve and endSolve
% 040102 hkh Return Iter = 1, n_r = n_J = 2
% 040602 med Help fix PriLev to PriLevOpt
% 041201 hkh weightType=3, dynamic weights, not possible, stop with error
% 041201 med Corrected typo
% 041202 hkh Use defblbu instead of most likely noncorrect code
% 041222 med Safeguard added for x_0
% 060818 med isnan checks removed for x_0
% 080607 hkh Change to iniSolveMini, endSolveMini
% 090717 med f calculations fixed
% 090910 hkh Add new output Iter, number of dual iterations in dlpdp-dwnnls
% 110722 hkh Remove check on weighting, weighting applied in llsAssign
