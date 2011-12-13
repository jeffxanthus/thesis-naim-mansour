% TOMLAB miqpBB MIQP/MILP, LP, QP Solver
%
% function Result = miqpBBTL(Prob)
%
% miqpBBTL converts the problem from the Tomlab structure format and
% calls either miqpBBs (sparse) or miqpBBd (dense).
% On return converts the result to the Tomlab structure format.
% Also see the help for miqpBB.m
%
% miqpBB solves the following mixed integer (linear or quadratic)
% programming problem (MILP, MIQP, LP, QP):
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
%   If F is empty, an LP or MILP problem is solved
%   Some or all x may be integer valued as specified by other input variables
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% Fields used in input structure Prob
% Use lpAssign, qpAssign or mipAssign to define the Prob structure,
%
% x_L, x_U  Lower and upper bounds on variables.
% b_L, b_U  Lower and upper bounds on linear constraints.
% A         Linear constraint matrix, dense or sparse m x n matrix.
% QP.c      Linear coefficients in objective function, size n x 1.
% QP.F      Quadratic matrix of size n x n.
% PriLevOpt Print Level.
%           (0 = off, 1 = summary, 2 = scalar information, 3 = verbose)
%           > 10 Pause statements, and maximal printing (debug mode).
%
% LargeScale If TRUE (=1), use sparse version, otherwise dense.
%
% MaxCPU    Maximum amount of time for solving, in seconds.
%
% -----------------------------------------------
% Fields used in Prob.optParam: (Structure with optimization parameters)
% -----------------------------------------------
% MaxIter   Limit of iterations
%
% -----------------------------------------------
% Prob.MIP      Structure with fields defining the integer properties of
%               the problem. The following fields are used:
%
%   IntVars:
%               If empty, all variables are assumed non-integer
%               If islogical(IntVars) (=all elements are 0/1), then
%               1 = integer variable, 0 = continuous variable.
%               If any element >1, IntVars is the indices for integer variables
%
%   VarWeight   Variable priorities. Lower value means higher priority.
%
%   fIP         An upper bound on the IP value wanted. Makes it possible to
%               cut branches and avoid node computations.
%
% -----------------------------------------------
% Fields used in Prob.DUNDEE:
% -----------------------------------------------
% kmax      Max dimension of reduced space (k), default n, set as 0 if LP.
% mlp       Maximum number of levels of recursion.
% stackmax  Maximum size of the LIFO stack storing info about B&B tree.
%           Default 5000.
%
% PrintFile Name of print file. Amount/print type determined by optPar(1).
%           Default name miqpBBout.txt.
%
% optPar    Vector of optimization parameters. If -999, set to default.
%           Length from 0 to 20 allowed.
%
% optPar(1) Print level in miqpBB:
%           = 0  Silent
%           = 1  Warnings and Errors
%           = 2  Summary information
%           = 3  More detailed information
%
% optPar(2)  tol   1E-10 Relative accuracy in solution
%            (Prob.optParam.eps_x)
% optPar(3)  emin  1.0   1.0 Use cscale (constraint scaling) 0.0 no scaling
% optPar(4)  sgnf  5E-4  Max rel error in two numbers equal in exact arithmetic
% optPar(5)  nrep  2     Max number of refinement steps
% optPar(6)  npiv  3     No repeat if no more than npiv steps were taken
% optPar(7)  nres  2     Max number of restarts if unsuccessful
% optPar(8)  nfreq 500   The max interval between refactorizations
%
%
% optPar(12) epsilon     Tolerance used for x value tests. DEF 1E-5
%                        Prob.optParam.eps_x used if present
%
% optPar(13) MIopttol    Tolerance used for function value tests 1E-4
%                        Prob.optParam.eps_f used if present
%
% optPar(14) fIP  Upper bound on f(x). Only consider solutions < fIP - MIopttol
%                 Prob.MIP.fIP used if present. DEFAULT infty = optPar(18) = 1E20
%
% optPar(15) timing   : If 1 - use timing, if 0 no timing (default)
% optPar(16) max_time : Maximal time allowed for the run in seconds
%            Default 4E3, i.e. 66 minutes
%
% optPar(17) branchType:  Branch on variable with highest priority. If tie:
%
%            = 1. Variable with largest fractional part, among those branch
%            on the variable giving the largest increase in the objective
%
%            = 2. Tactical Fletcher (PMO) branching. The var that solves
%            max(min(e+,e-)) is chosen. The problem than corresponding to
%            min(e+,e-) is placed on the stack first.
%
%            = 3. Tactical branching, Padberg/Rinaldi,91, Barahona et al.,89
%            (i) Choose the branching variable the one that most violates the
%            integrality restrictions. i.e. find  max(i){min(pi+,pi-)}
%            pi+ = int(x(i)+1) - x(i) , pi- = x(i) - int(x(i))
%            (ii) among those branch on the variable that gives the greatest
%            increase in the obj. function (iii) Finally a LOWER BOUND is
%            computed on the branched problems using the bounding method of
%            Fletcher and Leyffer (dual active set step) DEFAULT = 1
%
% optPar(18) ifsFirst: If 1, then only search for first ifs (ifail=6),DEFAULT 0
%
% optPar(19) infty       Real value for infinity  (default 1E20)
%
% Further fields used in Prob.DUNDEE:
%
% morereal   Number of extra REAL workspace locations. Set to <0
%            for problem dependent default value.
%
% moreint    Number of extra INTEGER workspace locations. Set to <0
%            for problem dependent default value.
%
% ----------------------------------------------------------------
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
%
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector (NOT USED).
% g_k      Exact gradient computed at optimum.
%
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
%
% ExitFlag - exit status from miqpBB.m (similar to TOMLAB)
% Inform   miqpBB information parameter
%           0 - Optimal solution obtained
%           1 - Error in parameters for BQPD
%           2 - Unbounded QP encountered
%           3 - Stack overflow NO ifs found
%           4 - Stack overflow some ifs obtained
%           5 - Integer infeasible
%           6 - (on I/O) only search for first ifs and stop
%           7 - Infeasible root problem
%
% rc       Reduced costs. NOT SET.
%
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
%
% Solver   Name of the solver (miqpBB).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.DUNDEE sub field
%
% kmax      Max dimension of reduced space (k), default n, set as 0 if LP.
% mlp       Maximum number of levels of recursion.
% stackmax  Maximum size of the LIFO stack storing info about B&B tree.
% mode      Mode of operation, default set as 2*Prob.WarmStart.
% x         Solution (Warm Start).
% k         Dimension of the reduced space (Warm Start).
% e         Steepest-edge normalization coefficients (Warm Start).
% ls        Indices of active constraints, first n-k used for warm start.
% lp        List of pointers to recursion information in ls (Warm Start).
% peq       Pointer to the end of equality constraint indices in ls (Warm Start).
% -----------------------------------------------
% Output fields in Result.MIP: (HKH: NOT USED NOW)
% -----------------------------------------------
% MIP.slack     Slack variables (m x 1 vector).
% MIP.ninf      Number of infeasibilities.
% MIP.sinf      Sum of infeasibilities.
% MIP.lpiter    Number of LP iterations.
% MIP.glnodes   Number of nodes visited.
% MIP.basis     Basis status of constraints + variables, (m + n x 1 vector).
%
% -----------------------------------------------------------------------
%
% For a problem description, see miqpAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written May 12, 2002.  Last modified Jun 7, 2008.

function Result = miqpBBTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1, error('miqpBBTL needs the Prob structure as input');end

global MAX_x % Max number of variables/constraints/resids to print

Prob.solvType = 11; % MIQP solver

Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='miqpBB';

PriLev = DefPar(Prob,'PriLevOpt',0);

% optPar vector. Fix length to 20 if smaller
optPar = DefPar(Prob.DUNDEE,'optPar',-999*ones(20,1));
if(length(optPar)<20), optPar(end+1:20) = -999; end

% Fill in some elements with data from other places if defaults are
% not given

eps_f = DefPar(Prob.optParam,'eps_f',1E-4);
eps_x = DefPar(Prob.optParam,'eps_x',1E-6);
fIP   = DefPar(Prob.MIP,'fIP',-999);

% Substitute with other values if defaults are specified in optPar
if(optPar(11)<0),     optPar(11) = eps_x; end
if(optPar(12)<0),     optPar(12) = eps_f; end
if(optPar(13)==-999), optPar(13) = fIP;   end

if ~isfield(Prob.optParam,'MaxIter')
   Prob.optParam.MaxIter=20000;
end

MaxCPU = DefPar(Prob,'MaxCPU',Inf);

if(isfinite(MaxCPU) & optPar(16)==-999 & optPar(15)~=0)
  % Enable timing if the user has not explicitly disabled it.
  optPar(16) = MaxCPU;
  optPar(15) = 1;
end

% Define lower and upper bound arrays for miqpBB
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints

% Set infty (BIG) = Real value for infinity, to default 1E20
if(optPar(19)<0), optPar(19) = 1E20; end

% Safe guard BIG, always >= 1E20
optPar(19) = max(1E20,optPar(19));
BIG        = optPar(19);

% Get single bounds vectors
% defblbu Order 0 gives [ bounds nonlinear linear]

[bl, bu, n, m] = defblbu(Prob, BIG);

% Maximum number of levels of recursion, typically 20, maximum = m+1.
mlp  = max(3,min(m+1,DefPar(Prob.DUNDEE,'mlp',40)));

% Maximum value of dimension of reduced space (k), default n, set as 0 if LP
% The check for LP/MILP is done later (Hzero)
kmax = min(n,DefPar(Prob.DUNDEE,'kmax',n ));

% Maximum value of LIFO stack (must be at least 100, check this)
stackmax = max(100,DefPar(Prob.DUNDEE,'stackmax',5000));

% Printfile
PrintFile = DefPar(Prob.DUNDEE,'PrintFile',[]);

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('miqpBBTL: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

if isempty(IntVars)
   MIP      = 0;
   Priority = [];
else
   MIP     = 1;
   VW      = DefPar(Prob.MIP,'VarWeight',[]);

   if ~isempty(VW)
      [i1,i2] = sort(VW);
      Priority(i2(i2)) = -i2;
      Priority = Priority - min(Priority);
   else
      Priority = [];
   end
end

if isempty(Prob.QP.F)
   Hzero = 1;
elseif n > 100
   Hzero = 0;
else
   % Only check for zero matrix on small problems
   Hzero = all(Prob.QP.F(:) == 0);
end

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m, error('Linear constraints A MUST have m rows!'); end 
end 

% Check if any linear part
c = Prob.QP.c(:);

% if isempty(c), c=zeros(n,1); end

Result.f_0=[];
% Determine type of problem
if isempty(c) | all(c==0)
   if Hzero
      kmax = 0; % 0 if LP problem
   end
   % Set up the constraint matrix A, transposed
   A = sparse([zeros(n,1),Prob.A']);
else
   if Hzero
      kmax = 0; % 0 if LP problem
   end
   % Set up the constraint matrix A, transposed 
   A = sparse([c(:),Prob.A']);
end


%if isempty(Prob.Name)
%   sprintf(ProbName,'Problem %d',Prob.P);
%else
%   ProbName=Prob.Name;
%end

% fLowBnd is fmin

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

moremem(1) = DefPar(Prob.DUNDEE,'morereal',-1);
moremem(2) = DefPar(Prob.DUNDEE,'moreint',-1);

if Prob.LargeScale
   [Inform, x_k, Obj, Iter] = miqpBBs(A, bl, bu, IntVars, Priority, ...
         'HxFunc', mlp, kmax, stackmax, optPar, PriLev, PrintFile, Prob,...
         moremem);
   
   Result.SolverAlgorithm='Sparse miqpBB MIQP/MILP code';
else
  [Inform, x_k, Obj, Iter] = miqpBBd(full(A), bl, bu, IntVars, Priority, ...
         'HxFunc', mlp, kmax, stackmax, optPar, PriLev, PrintFile, Prob,...
	 moremem);

   Result.SolverAlgorithm='Dense miqpBB MIQP/MILP code';
end

% Inform output ...
% f is either f or constraint deviation
% g is gradient is solution exists

switch Inform
   case {3,4,99}
     ExitFlag=1;  % Too many iterations or time limit
   case 2
     ExitFlag=2;  % Unbounded
   case {5,7}
     ExitFlag=4;  % Infeasible
   %case {x}
   %  ExitFlag=3;  % Rank problem
   case {1}
     ExitFlag=10; % Input errors
   otherwise
     % both 0 and 6 OK. 6 is search only for 1st ifs
     ExitFlag=0;
end

Result.f_k = Obj;
Result.x_k = x_k;

if ~isempty(c)
   if ~Hzero
      Result.g_k = Prob.QP.F*x_k + c;
   else
      Result.g_k=c;
   end
elseif Hzero
   Result.g_k=[];
else
   Result.g_k = Prob.QP.F*x_k;
end

if ~Hzero
   Result.GradEv   = Result.Iter;
else
   Result.GradEv   = 0;
end

Result.FuncEv    = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

%Result.MIP.slack=slack;
%Result.MIP.ninf=ninf;
%Result.MIP.sinf=sinf;
%Result.MIP.lpiter=lpiter;
%Result.MIP.glnodes=glnodes;
%Result.MIP.basis=basis;

if mA > 0
   Ax = Prob.A * x_k;
else
   Ax = [];
end

% Compute Result.xState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [],Prob.optParam.xTol, ...
                  Prob.optParam.bTol, [], bl, bu, 1); 

switch Inform
   case 0
    if(MIP)
      Text = 'Optimal ifs obtained';
    else
      Text = 'Optimal solution obtained';
    end
   case 1
      Text = 'Error in parameters for BQPD';
   case 2
      Text = 'Unbounded QP encountered';
   case 3
      Text = 'Stack overflow NO ifs found';
   case 4
      Text = 'Stack overflow some ifs obtained';
   case 5
      Text = 'Integer infeasible';
   case 6
      Text = '(on I/O) only search for first ifs and stop';
   case 7
      Text = 'Infeasible root problem';
   case 99
      Text = 'Time limit exceeded';
   otherwise
      Text = 'NOTE: UNKNOWN miqpBB Inform value.';
      % >8 = possible use by later sparse matrix codes
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nmiqpBB solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('miqpBB: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   if MIP
      fprintf('Iterations   %7d. ',Iter);
   else
      fprintf('Nodes visited%7d. ',Iter);
   end
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
% 020621 hkh  New interface written
% 020810 hkh  Revision
% 030114 ango Add iprint in optPar(1)
% 030123 ango Correct comments about stackmax
% 030128 ango Revised comments
% 030129 ango Change name: miqpbb->miqpBB
% 030206 hkh  Expand optPar to 20
% 030220 hkh  Restrict mlp and kmax,  min(m, mlp), min(n,kmax), mlp default 40
% 030224 ango Fix optPar error
% 030225 ango Change min value for mlp to m+1 or 3
% 031113 hkh  Change infty,BIG to 1E12 instead of 1E20, set default optPar(19)
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040102 hkh  Return Iter, all others 0
% 040803 med  Added pragmas for MATLAB Compiler
% 041202 hkh  filterSQP fails if BIG < 1D20, correct this bug, safe guard BIG
% 041202 hkh  Revise calls to defblbu and StateDef
% 041210 ango Add support for Prob.MaxCPU
% 060306 med  v_k removed from help
% 070222 hkh  Revised IntVars handling, use new format
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
