% TOMLAB BQPD QP/LP Solver
%
% function Result = bqpdTL(Prob)
%
% BQPD solves the following linear or quadratic
% programming problem (LP, QP):
%
%   minimize   0.5 * x'*F*x + c'x     subject to:
%      x             x_L <=    x   <= x_U
%                    b_L <=   Ax   <= b_U
%   where
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the Dundee solvers sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense.
%   If F is empty, an LP problem is solved
%
%   If F is not positive semidefinite, only a local solution is found
%
% INPUT:
% Prob   Problem structure in Tomlab format
%
% Fields used in input structure Prob
% Use qpAssign to define the Prob structure,
% (or call Prob=ProbDef; and define each field)
%
% x_L, x_U   Bounds on variables.
% b_L, b_U   Bounds on linear constraints.
% A          Linear constraint matrix.
% QP.c       Linear coefficients in objective function.
% QP.F       Quadratic matrix of size n x n.
% PriLevOpt  Print Level.
%            (0 = off, 1 = summary, 2 = scalar information, 3 = verbose)
% WarmStart  If TRUE (=1), use warm start, otherwise cold start.
% LargeScale If TRUE (=1), use sparse version, otherwise dense.
%
% -----------------------------------------------
% Fields used in Prob.DUNDEE:
% -----------------------------------------------
%  QPmin    Lower bound for the QP subproblems.
%           If not set, use Prob.f_Low. Default: -1E300
%
% callback  If 1, use a callback to Matlab to compute QP.F * x for different x.
%           Faster when F is very large and almost dense, avoiding
%           copying of F from Matlab to MEX.
% kmax      Max dimension of reduced space (k), default n, set as 0 if LP.
% mlp       Maximum number of levels of recursion.
% mode      Mode of operation, default set as 2*Prob.WarmStart.
%
% x         Solution (Warm Start).
% k         Dimension of the reduced space (Warm Start).
% e         Steepest-edge normalization coefficients (Warm Start).
% ls        Indices of active constraints, first n-k used for warm start.
% lp        List of pointers to recursion information in ls (Warm Start).
% peq       Pointer to the end of equality constraint indices in ls (Warm Start).
% PrintFile Name of print file. Amount/print type determined by optPar(1).
%           Default name bqpd.txt
%
% optPar    Vector of optimization parameters. If -999, set to default.
%           Length from 0 to 20 allowed. Elements used:
%
% optPar(1): iprint 0     Print level in PrintFile
% optPar(2): tol    1E-10 Relative accuracy in solution
% optPar(3): emin   1.0   1.0 Use cscale (constraint scaling) 0.0 no scaling
% optPar(4): sgnf   5E-4  Max rel error in two numbers equal in exact arithmetic
% optPar(5): nrep   2     Max number of refinement steps
% optPar(6): npiv   3     No repeat if no more than npiv steps were taken
% optPar(7): nres   2     Max number of restarts if unsuccessful
% optPar(8): nfreq  500   The max interval between refactorizations
% optPar(19):infty  1E20  A large value representing infinity
%
%
% --------------------------------------------
%
% OUTPUT:
% Result     Structure with results (see ResultDef.m):
%
% f_k        Function value at optimum or constraint deviation if infeasible.
% x_k        Solution vector.
% x_0        Initial solution vector.
% g_k        Exact gradient computed at optimum.
%
% xState     State of variables. Free == 0; On lower == 1; On upper == 2;
%            Fixed == 3;
% bState     State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%            Equality == 3;
% v_k        Lagrangian multipliers (for bounds + dual solution vector).
%
% ExitFlag   Exit status from bqpd.m (similar to Tomlab).
% Inform     BQPD information parameter.
%            0 - Solution obtained
%            1 - Unbounded problem detected (f(x)<=fLow occurred)
%            2 - Lower bound bl(i) > bu(i) (upper bound) for some i
%            3 - Infeasible problem detected in Phase 1
%            4 - Incorrect setting of m, n, kmax, mlp, mode or tol
%            5 - Not enough space in lp
%            6 - Not enough space for reduced Hessian matrix (increase kmax)
%            7 - Not enough space for sparse factors
%            8 - Maximum number of unsuccessful restarts taken
%
% Iter       Number of iterations.
% MinorIter  Number of minor iterations. Always set to 0.
% FuncEv     Number of function evaluations. Set to Iter.
% GradEv     Number of gradient evaluations. Set to Iter.
% ConstrEv   Number of constraint evaluations. Set to 0.
%
% QP.B       Basis vector in Tomlab QP standard.

% Solver           Name of the solver (BQPD).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.DUNDEE sub field:
%
% kmax      Max dimension of reduced space (k), default n, set as 0 if LP.
% mlp       Maximum number of levels of recursion.
% mode      Mode of operation, default set as 2*Prob.WarmStart.
% x         Solution (Warm Start).
% k         Dimension of the reduced space (Warm Start).
% e         Steepest-edge normalization coefficients (Warm Start).
% ls        Indices of active constraints, first n-k used for warm start.
% lp        List of pointers to recursion information in ls (Warm Start).
% peq       Pointer to the end of equality constraint indices in ls (Warm Start).
%
% -----------------------------------------------------------------------
%
% For a problem description, see qpAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written May 12, 2002.  Last modified Jul 17, 2009.

function Result = bqpdTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1, error('bqpdTL needs the Prob structure as input');end

global MAX_x % Max number of variables/constraints/resids to print

Prob.solvType = 2; % QP solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='BQPD';

LargeScale = DefPar(Prob,'LargeScale',0);

switch(LargeScale)
case 1,
   Result.SolverAlgorithm = 'Sparse Active Set Method';
otherwise,
   Result.SolverAlgorithm = 'Dense Active Set Method';
end

PriLev = DefPar(Prob,'PriLevOpt',0);

optPar    = DefPar(Prob.DUNDEE,'optPar',-999*ones(20,1));
if(length(optPar)<20), optPar(end+1:20)=-999; end

% Set infty (BIG) = Real value for infinity, to default 1E20
if(optPar(19)<0), optPar(19) = 1E20; end

% Safe guard BIG, always >= 1E20
optPar(19) = max(1E20,optPar(19));
BIG        = optPar(19);

% Define lower and upper bound arrays for BQPD
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints

[bl, bu, n, m] = defblbu(Prob, BIG,1);

% Initial checks on the inputs

fLow = DefPar(Prob.DUNDEE,'QPmin',Prob.f_Low);

nb = n+m;

% Maximum number of levels of recursion, typically 20, maximum = m+1.
mlp = max(3,min(m+1,DefPar(Prob.DUNDEE,'mlp',40)));

% Maximum value of dimension of reduced space (k), default n, set as 0 if LP
kmax = min(n,DefPar(Prob.DUNDEE,'kmax',n ));

PrintFile = DefPar(Prob.DUNDEE,'PrintFile',[]);
callback  = DefPar(Prob.DUNDEE,'callback',0);

% Check on Warm start

if Prob.WarmStart
   % Warm start for DUNDEE solver
   x_0  = Prob.DUNDEE.x;
   e    = Prob.DUNDEE.e;
   ls   = Prob.DUNDEE.ls;
   lp   = Prob.DUNDEE.lp;
   peq  = Prob.DUNDEE.peq;
   k    = Prob.DUNDEE.k;
   if isfield(Prob.DUNDEE,'mode')
      mode = Prob.DUNDEE.mode;
   else
      mode = 2;
   end
else
   % Initial values
   x_0     = DefPar(Prob,'x_0',zeros(n,1));
   % Safe guard
   x_0     = max(bl(1:n),min(x_0,bu(1:n)));
   e    = [];
   ls   = [];
   lp   = [];
   peq  = [];
   k    = [];
   if isfield(Prob.DUNDEE,'mode')
      mode = Prob.DUNDEE.mode;
   else
      mode = 0;
   end
end

if Prob.WarmStart
   Hzero = 0;
else
   if isempty(Prob.QP.F)
      Hzero = 1;
   elseif n > 100
      Hzero = 0;
   else
      % Only check for zero matrix on small problems
      Hzero = all(Prob.QP.F(:) == 0);
   end
end

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   if mA~=m, error('Linear constraints A MUST have m rows!'); end
end

% Check if any linear part
c = Prob.QP.c(:);

% Determine type of problem
if isempty(c) | all(c==0)
   if Hzero
      Result.f_0 = 0;
      kmax = 0; % 0 if LP problem
   else
      Result.f_0=0.5*(x_0'*Prob.QP.F*x_0);
   end
   % Set up the constraint matrix A, transposed
   A = sparse([zeros(n,1),Prob.A']);
else
   if Hzero
      if Prob.WarmStart
         Result.f_0=c(1:n)'*x_0;
      else
         Result.f_0=c(1:n)'*x_0;
         %Result.f_0=0;
      end
      kmax = 0; % 0 if LP problem
   else
      Result.f_0=0.5*(x_0'*Prob.QP.F*x_0) + c(1:n)'*x_0;
   end
   % Set up the constraint matrix A, transposed
   A = sparse([c(:),Prob.A']);
end

if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Number of variables          %d.\n',n)
      fprintf('Number of linear constraints %d.\n',m)
   end
   if PriLev >= 1000
      xprinte(bl(1:n),'x_L: ');
      xprinte(bu(1:n),'x_U: ');
      xprinte(bl(n+1:nb),'b_L: ');
      xprinte(bu(n+1:nb),'b_U: ');
      xprinte(c,'c: ');
      if PriLev > 0, pause; end
   end
   if PriLev >= 7
      if PriLev >= 10, pause; end
      PrintMatrix(full(Prob.A),'Constraint matrix A:')
      if PriLev >= 10, pause; end

      disp('x_L x_0 x_U');
      if Prob.WarmStart
         mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      else
         mPrint([bl(1:n) bu(1:n)],'xl/xu:')
      end
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 10, pause; end
   end
end

%if isempty(Prob.Name)
%   sprintf(ProbName,'Problem %d',Prob.P);
%else
%   ProbName=Prob.Name;
%end

% info(1) = Number of pivots, i.e. number of iterations
% k         Dimension of the reduced space
% kmax      Maximum value of k (kmax = 0 if LP)
% x         Solution
% e         Steepest-edge normalization coefficients
% ls        Indices of active constraints, first n-k used for warm start
% mlp       Maximum number of levels of recursion
% peq       Pointer to the end of equality constraint indices in ls
% mode      Mode of operation
% Iter      info(1) - assume length of info is 1 !!!

% Warm Start if mode >= 2

% Create A from Prob.A and c, as in MINOS

% moremem
moremem(1) = DefPar(Prob.DUNDEE,'morereal',-1);
moremem(2) = DefPar(Prob.DUNDEE,'moreint',-1);

if callback
   if LargeScale
      [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpds(A, x_0, ...
         bl, bu, 'HxFunc',fLow,mlp,mode,kmax,PriLev,PrintFile, ...
         k,ls,e,peq,lp,optPar,Prob,moremem);
   else
      [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpdd(full(A),...
         x_0,  bl, bu, 'HxFunc',fLow,mlp,mode,kmax,PriLev,PrintFile, ...
         k,ls,e,peq,lp,optPar,Prob,moremem);
   end

elseif LargeScale

   [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpds(A, x_0, ...
      bl, bu, Prob.QP.F,fLow,mlp,mode,kmax,PriLev,PrintFile, ...
      k,ls,e,peq,lp,optPar,Prob,moremem);
else

   [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpdd(full(A),...
      x_0,bl, bu, Prob.QP.F,fLow,mlp,mode,kmax,PriLev,PrintFile, ...
      k,ls,e,peq,lp,optPar,Prob,moremem);
end

% f is either f or constraint deviation

Result.DUNDEE.x     = x_k;

% Warm start
Result.DUNDEE.k     = k;
Result.DUNDEE.ls    = ls;
Result.DUNDEE.e     = e;
Result.DUNDEE.peq   = peq;
Result.DUNDEE.lp    = lp;

% e         Steepest-edge normalization coefficients
% ls        Indices of active constraints, first n-k used for warm start
% mlp       Maximum number of levels of recursion
% peq       Pointer to the end of equality constraint indices in ls

switch Inform
   case 8
     ExitFlag=1;  % Too many iterations
   case 1
     ExitFlag=2;  % Unbounded
   case 3
     ExitFlag=4;  % Infeasible
   %case {x}
   %  ExitFlag=3;  % Rank problem
   case {2,4,5,6,7}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = v_k;


% Gradient
if kmax == 0
   Result.g_k=c;
else
   Result.g_k=g;
end

Result.FuncEv    = 0;
Result.GradEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

if mA > 0
   Result.Ax = Prob.A * x_k;
   Result = StateDef(Result, x_k, Result.Ax, [],Prob.optParam.xTol, ...
                     Prob.optParam.bTol, [], bl, bu, 1);
else
   % Compute Result.xState and Result.QP.B only
   Result = StateDef(Result, x_k, [], [],Prob.optParam.xTol, ...
                     [], [], bl, bu, 1);
end


switch Inform
   case 0
      Text = 'Solution obtained';
   case 1
      Text = 'Unbounded problem detected (f(x)<=fLow occurred)';
   case 2
      Text = 'Lower bound bl(i) > bu(i) (upper bound) for some i';
   case 3
      Text = 'Infeasible problem detected in Phase 1';
   case 4
      Text = 'Incorrect setting of m, n, kmax, mlp, mode or tol';
   case 5
      Text = 'Not enough space in lp';
   case 6
      Text = 'Not enough space for reduced Hessian matrix (increase kmax)';
   case 7
      Text = 'Not enough space for sparse factors';
   case 8
      Text = 'Maximum number of unsuccessful restarts taken';
   otherwise
      Text = 'NOTE: UNKNOWN BQPD Inform value.';
      % >8 = possible use by later sparse matrix codes
end

Result.ExitText = Text;


if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nBQPD solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('BQPD: Inform = %2d, ',Inform)
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
         MAX_x=length(x_k);
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(length(x_k),MAX_x)),'x_k:  ');
   end

end
Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 020512 hkh  New interface written
% 020525 hkh  Further modifications
% 020618 hkh  Revision of comments, efficient use, final parameters
% 020630 hkh  Using callback function HxFunc. Use Prob.x_0 for initial x
% 020701 hkh  Compute Result.f_0 correctly for cold start
% 030106 hkh  Different SolverAlgorithm text dependent on sparse/dense version
% 030113 hkh  New way to set fLow, correct comments
% 030128 ango Revised comments
% 030129 ango Safe guard x_0
% 030206 hkh  Increase optPar to 20
% 030220 hkh  Restrict mlp and kmax,  min(m, mlp), min(n,kmax)
% 030225 ango Change min. value for mlp to m+1 or 3
% 031114 hkh  Change infty,BIG to 1E12 instead of 1E20, set default optPar(19)
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040102 hkh  Return Iter, all others 0
% 040607 med  Removed double assign for x_0
% 040607 med  f_0 print out removed
% 040803 med  Added pragmas for MATLAB Compiler
% 041202 hkh  Revise calls to defbl/StateDef, return Result.Ax
% 060818 hkh  Use Prob.f_Low instead of -1E300 defining fLow, if QPmin not set
% 060818 med  isnan checks removed for x_0
% 080607 hkh  Use iniSolveMini and endSolveMini
