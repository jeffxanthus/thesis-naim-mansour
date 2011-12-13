% TOMLAB TLSQR Least Squares Solver
%
% LSQR  finds a solution x to the following problems:
%
% 1. Unsymmetric equations - solve  A*x = b
%
% 2. Linear least squares  - solve  A*x = b in the least-squares sense
%
% 3. Damped least squares  - solve  (   A    )*x = ( b )
%                                   ( damp*I )     ( 0 ) in least-squares sense
%
% where A is a matrix with m rows and n columns, b is an m-vector, and damp is
% a nonnegative scalar (All quantities are real). The matrix A is intended to
% be a large and sparse array, but both full and sparse Matlab arrays are
% handled.
%
% Do   help Tlsqr  to see a discussion of the numerical aspects of using the
%                  LSQR algorithm implemented in the solver Tlsqr
%
% Note! LSQR is not intended to solve problems with linear constraints.
% If defined, LSQR tries to add the linear constraints as additional residuals
% but an optimal result with low accuracy is only obtained if all linear
% constraints are equality constraints.
% To solve large-scale linear least squares with linear constraints, the solver
% pdco can be used. It calls Tlsqr as a subsolver.
%
% function Result = TlsqrTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints (avoid giving these, see Note).
% A         Linear constraint matrix (avoid giving these, see Note).
% LS.C      Linear matrix m x n.
% LS.y      Data vector m x 1.
% LS.damp   Damping parameter >= 0. NOTE!!! If damp < 0 set to abs(damp).
% LS.condLim An upper limit on cond(Abar), where Abar = [A;damp * I].
% PriLevOpt Print Level.
% optParam.MaxIter - MaxIter in Tlsqr
% optParam.xTol    - aTol in Tlsqr
% optParam.bTol    - bTol in Tlsqr
%
% SOL.PrintFile  Name of print file. For example:
%                Prob.SOL.PrintFile = 'tlP.txt'
%
% x_0       Initial solution vector. Normally LSQR starts with x=0, but
%           if x_0 is given (n-vector), LSQR tries a warm start with x_0.
%
% D         Preconditioning. Only n diagonal elements of D are given.
%
% Alloc     Advanced memory handling if Alloc > 0 (default empty (or 0))
%           See help Tlsqr.m , input parameter m, for a description of
%           how to use Alloc.
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% r_k      Residual vector.
% J_k      Jacobian, is just the Prob.LS.C matrix.
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial  solution vector, empty if not given as a n-vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status  (similar to TOMLAB).
% Inform   LSQR information parameter.
%
%          0 = Both equality and inequality constraints are compatible
%              and have been satisfied. or
%              x = 0  is the exact solution.  No iterations were performed.
%          1 = The equations A*x = b are probably compatible.
%              Norm(A*x - b) is sufficiently small
%              given the values of aTol and bTol.
%          2 = damp is zero. The system A*x = b is probably not compatible.
%              A least-squares solution has been obtained that is
%              sufficiently accurate,  given the value of aTol.
%          3 = damp is nonzero. A damped least-squares solution has been
%              obtained that is sufficiently accurate, given the value of aTol
%          4 = An estimate of cond(Abar) has exceeded condLim. The system
%              A*x = b appears to be ill-conditioned.  Otherwise, there
%              could be an error in the internal matrix product routine.
%          5 = The iteration limit MaxIter was reached.
%
% Iter     Number of iterations, set to -1.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% Solver   Name of the solver (Tlsqr).
% SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------------------------------
%
% For a problem description, see Tlsqr.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1994-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Nov 8, 2000.    Last modified Jun 6, 2008.

function Result = TlsqrTL(Prob)

if nargin < 1, error('TlsqrTL needs the Prob structure as input');end

global MAX_x % Max number of variables/constraints/resids to print

Prob.solvType = 5; % LLS solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='Tlsqr';
Result.SolverAlgorithm='TOMLAB LSQR';

PriLev=Prob.PriLevOpt;

% Initial checks on the inputs
[bl,bu,n,mA] = defblbu(Prob,inf,1);

% Linear least squares
C = Prob.LS.C;
y = Prob.LS.y(:);
m = length(y);

x_0     = Prob.x_0(:);
if length(x_0) < n, x_0=[]; end

% Safe-guard starting point
x_0   = max(bl(1:n),min(x_0,bu(1:n)));

x_L = bl(1:n);
x_U = bu(1:n);

if mA > 0
   nA = size(Prob.A,2);
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   b_L = bl(n+1:n+mA);
   b_U = bu(n+1:n+mA);
else
   b_L = [];
   b_U = [];
end

% Set up the constraint matrices E x = f and G x >= h
% Weight for equalities
w   = 1000;
wLU = 100;

if mA > 0
   ix1 = find( b_L == b_U & ~isinf(b_L) );
   ix2 = find( b_L ~= b_U & ~isinf(b_U) );
   ix3 = find( b_L ~= b_U & ~isinf(b_L) );
   mE  = length(ix1);
   if mE == 0
      E = zeros(0,n);
   else
      E   = w*Prob.A(ix1,:);
   end
   f   = w*b_L(ix1);
   h   = w*[-b_U(ix2); b_L(ix3)];
   mG  = length(h);
   if mG > 0
      % Create G, with slack
      G   = w*[[-Prob.A(ix2,:); Prob.A(ix3,:)],speye(mG,mG)];
   else
      G   = zeros(0,n);
   end
else
   E   = zeros(0,n);
   f   = [];
   G   = zeros(0,n);
   h   = [];
   mE  = 0;
   mG  = 0;
end
ixL = find(~isinf(x_L));
ixU = find(~isinf(x_U));

mL = length(ixL);
if mL > 0
   GL = sparse(1:mL,ixL,wLU*ones(mL,1),mL,n+mG);
   G  = [G,sparse(mG,mL);GL,wLU*speye(mL,mL)];
   h  = [h;wLU*x_L(ixL)];
end
mU = length(ixU);
if mU > 0
   GU = sparse(1:mU,ixU,-wLU*ones(mU,1),mU,n+mG+mL);
   G  = [G,sparse(mG+mL,mU);GU,wLU*speye(mU,mU)];
   h  = [h;-wLU*x_U(ixU)];
end
k = size(G,2)-n;

zE = 100;
zG = 100;
if issparse(C)
   A = sparse([C,sparse(m,k);zE*E,sparse(mE,k);zG*G]);
else
   A = full([C,zeros(m,k);zE*E,zeros(mE,k);zG*G]);
end
b = [y;zE*f;zG*h];

% Linear least squares
if isempty(x_0)
   r = -y;
   Result.f_0=0.5*(r'*r);
else
   r = C*x_0-y;
   Result.f_0=0.5*(r'*r);
end

% These parameters should be sent some way
if isfield(Prob.LS,'damp')
   damp    = abs(Prob.LS.damp);
else
   damp    = 0;
end
aTol       = Prob.optParam.xTol;
bTol       = Prob.optParam.bTol;
if isfield(Prob.LS,'condLim')
   condLim = Prob.LS.condLim;
else
   condLim = [];
end
%MaxIter    = 100*(n+k);
MaxIter    = Prob.optParam.MaxIter;
WantStdErr = 1;
D          = [];

% ANGO - to get correct value of MwAlloc, I needed to switch
% order between the following row...
[M,N]=size(A);

% --- and this block of code ---
Alloc      = DefPar(Prob,'Alloc',[]);
if ~isempty(Alloc)
   MwAlloc = [0;Alloc];
else
   MwAlloc = M;
end
% --- end here ----


if ~isempty(x_0)
   if length(x_0) < N
      x_0 = [x_0;zeros(N-length(x_0),1)];
   end
end

% Can set nOut to a filename, or Fortran file unit (writes to lsqrout.txt)
nOut = DefPar(Prob.SOL,'PrintFile',[]);

[ x_k, Inform, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
    Tlsqr( MwAlloc, N, A, [], [], b, damp, aTol, bTol, condLim, ...
    MaxIter, WantStdErr, nOut, D, x_0 );

s = x_k(n+1:n+k);
x_k = x_k(1:n);

if mA > 0 & Inform <= 3
   Ax = Prob.A*x_k;
else
   Ax = [];
end

% Get the Lagrange multipliers back to the correct place in TOMLAB
%if length(v) > 0
%   cLamda = zeros(n+m,1);
%   CLamda(1:n) = v(1:n);
%   cLamda(n+ix1) = v(n+1:n+length(ix1));
%   j = n+length(ix1)+length(ix2);
%   cLamda(n+ix2) = v(n+length(ix1)+1:j);
%
%   for i = 1:length(ix3)
%       if v(j+i) ~= 0
%          cLamda(n+ix3(i)) = v(j+i);
%       end
%   end
%end

switch Inform
   case {0,1,2,3}
     ExitFlag=0;  % Convergence
   case 5
     ExitFlag=1;  % Too many iterations
   case {4}
     ExitFlag=3;  % Rank problem
%  case {5}
%    ExitFlag=10; % Input errors
%  case -1
%    ExitFlag=2;  % Unbounded
%  otherwise
%    ExitFlag=4;  % Infeasible
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

Result.FuncEv    = 0;
Result.ResEv     = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState, Result.bState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [], Prob.optParam.xTol, ...
         Prob.optParam.bTol, [], bl, bu, 1);

switch Inform
   case 0
      Text = 'x = 0  is the exact solution.  No iterations were performed.';
   case 1
      Text = str2mat(...
           'The equations A*x = b are probably compatible.', ...
           'Norm(A*x - b) is sufficiently small', ...
           'given the values of aTol and bTol.');
   case 2
      Text = str2mat( ...
           'damp is zero. The system A*x = b is probably not compatible.',...
           'A least-squares solution has been obtained that is', ...
           'sufficiently accurate,  given the value of aTol.');
   case 3
      Text = str2mat( ...
        'damp is nonzero. A damped least-squares solution has been ',...
        'obtained that is sufficiently accurate, given the value of aTol');
   case 4
      Text = str2mat( ...
        'An estimate of cond(Abar) has exceeded condLim. The system', ...
        'A*x = b appears to be ill-conditioned.  Otherwise, there', ...
        'could be an error in the internal matrix product routine.');
   case 5
      Text = 'The iteration limit MaxIter was reached.';
   otherwise
      Text = 'Tlsqr: System error, illegal return parameter';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLSQR solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LSQR: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   fprintf('Sum of squares at solution x %26.18f\n\n',Result.f_k);
   fprintf('Iterations %11d\n',Iter);
   fprintf('Variables  %11d\n',n);
   fprintf('Slacks     %11d\n',k);
   fprintf('Equations  %11d\n',size(A,1));
   fprintf('rNorm   %13.7e\n',rNorm);
   fprintf('xNorm   %13.7e\n',xNorm);
   fprintf('aNorm   %13.7e\n',aNorm);
   fprintf('aCond   %13.7e\n',aCond);
   fprintf('arNorm  %13.7e\n',arNorm);
   fprintf('\n');

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=length(x_k);
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:   ');
   end
   if PriLev > 1 & WantStdErr
      fprintf('Standard Errors for x = \n');
      xprinte(StdErr(1:min(n,MAX_x)),'Std: ');
   end
   if PriLev > 2
      fprintf('Slacks added for the inequalities and bounds \n');
      fprintf('( > 0 if OK)  = \n');
      xprinte(s(1:min(k,MAX_x)),'s:   ');
   end
end
Result=endSolveMini(Prob,Result);

% MODIFICATION LOG
%
% 030119  hkh  Revision for Tomlab v4.0, change name and call of mex
% 030124  hkh  Use Prob.Alloc, if defined
% 030221  hkh  x should be x_k in length(x_k). nOut was undefined
% 030918  ango Fix misspelling of MwALLOC (should be MwAlloc)
% 040102  hkh  Revision for v4.2, call iniSolve and endSolve
% 040102  hkh  Return Iter , n_r = n_J = Iter
% 040429  hkh  Safe guard y in LS.y to be a column vector
% 040429  hkh  Use intermediate sparse and speye instead of zeros and eye
% 040602  med  Help fix PriLev to PriLevOpt
% 041123  hkh  Must set NwAlloc to M, not m
% 041203  hkh  Revise calls to defblbu and StateDef, use output from defblbu
% 041203  hkh  Suspect bugs in previous code when x_L, x_U had infs
% 041222  med  Safeguard added for x_0
% 050801  med  return removed from input check
% 051216  med  Help updated with Inform values and text
% 051216  med  Additional text removed
% 060730  hkh  Cleaning up. Use abs(damp) to ensure >= 0
% 060818  med  isnan checks removed for x_0
% 070529  ango Change nOut setting, now uses Prob.SOL.PrintFile if set
% 080606  med  Switched to iniSolveMini
% 080607  hkh  Switched to endSolveMini
