% TOMLAB PENSDP SDP Solver
%
% function Result = pensdpTL(Prob)
%
% -----------------------------------------------------------------------------
% INPUT:
% Prob   Problem structure in Tomlab format
%
%  x_L, x_U   Bounds on variables.
%  b_L, b_U   Bounds on linear constraints.
%  A          Linear constraint matrix.
%  QP.c       Linear coefficients in objective function.
%
%  PriLevOpt  Print Level in pensdpTL and MEX interface
%             (0 = off, 1 = summary, 2 = scalar information, 3 = verbose).
%
% -----------------------------------------------
% Fields used in Prob.PENOPT:
% -----------------------------------------------
%
%  SDP          Structure array with Linear Matrix Inequality (LMI) data.
%               Please refer to Tomlab User's Guide for information on
%               how to set this field.
%
%
%  ioptions     (15x1)-vector with integer options.
%               Set any element < 0 (except ioptions(1)) to make
%               standard Tomlab parameter take precedence.
%               Where applicable, standard Tomlab parameter is given in [].
%               Default values are PENSDP defaults and are used
%               if nothing else is provided.
%
%  ioptions(1)  Set to 0 for PENSDP specified default values
%               in ioptions and foptions.
%
%  ioptions(2)  Maximum number of iterations of the overall algorithm.
%               [Prob.optParam.MaxIter] Default: 50
%
%  ioptions(3)  Maximum number of iterations for unconstrained minimization.
%               [Prob.optParam.MinorIter] Default: 100
%
%  ioptions(4)  Print level in PENSDP solver: 0=silent, 1,2,3=summary,brief,full
%               Default: 0
%
%  ioptions(5)  Hessian density check: 0=automatic check, 1=assume dense,
%               2=assume sparse
%               Default: 0 (no)
%
%  ioptions(6)  Use/do not use linesearch in unconstrained subproblem
%               Default: 0 (no)
%
%  ioptions(7)  Write solution vector to output file
%               Default: 0 (no)
%
%  ioptions(8)  Write computed multipliers to output file
%               Default: 0 (no)
%
%  ioptions(9)  Mode of solution of the Newton system.
%                0  Cholesky method (default)
%                1  Preconditioned conjugate gradient method
%                2  Preconditioned CG method with approximate
%                   Hessian calculation
%                3  Preconditioned CG method with exact Hessian not explicitly
%                   assembled
%                4  Symmetric Gauss-Seidel
%
%  ioptions(10) Preconditioner type for the CG method.
%                0  No preconditioner (default)
%                1  Diagonal
%                2  LBFGS
%                3  Appoximate inverse
%                4  Symmetric Gauss-Seidel
%
%  ioptions(11) Print DIMACS error measures.
%                Default: 0 (no)
%
%  ioptions(12) Trust-Region mode.
%                0  Modified Newton method (default)
%                1  Trust Region method
%
%  ioptions(13) Hybrid Mode - switch to Cholesky if PCG causes trouble.
%                0  No (default)
%                1  Yes, use diag. preconditioner after Cholesky step
%                2  Yes, use inverse Cholesky prec. after Cholesky step
%
%  ioptions(14) Auto Init
%                0 No
%                1 Yes (Default)
%
%  ioptions(15) Stop Mode
%                0 Strict
%                1 Heuristic (default)
%
%
%  foptions     (12x1)-vector, floating point parameters
%
%  foptions(1)  Scaling factor for linear constraints
%               Must be positive. Default: 1.0
%
%  foptions(2)  Restriction for multiplier update for linear constraints
%               Default: 0.7
%
%  foptions(3)  Restriction for multiplier update for matrix constraints
%               Default 0.1
%
%  foptions(4)  Stopping criterium for overall algorithm (rel. change in f)
%               [Prob.optParam.eps_f] Default: 1.0e-7
%
%  foptions(5)  Lower bound for penalty parameters
%               Default: 1.0e-6
%
%  foptions(6)  Lower bound for multipliers0
%               Default: 1.0e-14
%
%  foptions(7)  Stopping criterium for unconstrained minimization
%               Default: 1.0e-2
%
%  foptions(8)  Initial penalty value. Set lower (0.01-0.1) to maintain
%               feasibility when starting from a feasible point.
%               Default: 1.0
%
%  foptions(9)  Penalty update; when set to 0.0, it is computed
%               automatically.
%               Default: 0.0
%
%  foptions(10) Update of alpha; should either be equal to 1.0
%               (= no update) or smaller than 1.0 (e.g. 0.7).
%               Default: 1.0
%
%  foptions(11) Precision of the KKT conditions.
%               Default: 1.0e-7
%
%  foptions(12) Stopping criterion for the conjugate gradient
%               algorithm.
%               Default: 5.0e-2
%
% -----------------------------------------------------------------------------
% OUTPUT:
% Result  Structure with results (see ResultDef.m):
%
%   f_k       Function value at optimum or constraint deviation if infeasible.
%   x_k       Solution vector.
%   x_0       Initial solution vector.
%   g_k       Exact gradient computed at optimum.
%
%   xState    State of variables. Free == 0; On lower == 1; On upper == 2;
%             Fixed == 3;
%   bState    State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%             Equality == 3;
%   v_k       Lagrangian multipliers (for bounds + dual solution vector).
%
%   ExitFlag  Exit status from pensdp.m (similar to Tomlab).
%   Inform    PENSDP information parameter.
%              0 - Solution obtained
%              1 - No progress in objective value, problem may be infeasible
%              2 - Cholesky factorization of Hessian failed. The result may still be useful
%              3 - Maximum iteration limit exceeded. The result may still be useful
%              4 - Linesearch failed. The result may still be useful
%              5 - Wrong input parameters (ioptions, foptions)
%              6 - Memory error
%              7 - Unknown error. Contact support@tomopt.com
%
%   Iter      Number of iterations.
%   FuncEv    Number of function evaluations. Set to Iter.
%   GradEv    Number of gradient evaluations. Set to Iter.
%   ConstrEv  Number of constraint evaluations. Set to 0.
%   QP.B      Basis vector in Tomlab QP standard.
%   MinorIter Number of minor iterations. Always set to 0.
%
%   Solver    Name of the solver (PENSDP).
%   SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------
% Fields returned in Result.PENOPT:
%
%  iresults     Integer results:
%
%   iresults(1)  Number of outer iterations
%   iresults(2)  Number of Newton steps
%   iresults(3)  Number of linesearch steps
%   iresults(4)  Elapsed time in seconds
%
%  fresults     Floating point results
%
%   fresults(1)  Relative precision at x_k
%   fresults(2)  Feasibility of linear inequalities at x_k
%   fresults(3)  Feasibility of matrix inequalities at x_k
%   fresults(4)  Complementary slackness of linear inequalities at x_k
%   fresults(5)  Complementary slackness of matrix inequalities at x_k
%
%  info      PENSDP info value
%
% -----------------------------------------------------------------------
%
% For a problem description, see sdpAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 3, 2002.  Last modified Jun 7, 2008.

function Result = pensdpTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1, error('pensdpTL needs the Prob structure as input');end

global MAX_x% Max number of variables/constraints/resids to print

Prob.solvType = 13; % LMI solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='PENSDP';
Result.SolverAlgorithm='LMI Solver PENSDP 2.2';

PriLev=DefPar(Prob,'PriLevOpt',0);

x_0     = Prob.x_0;

if isfield(Prob.PENOPT,'pen')
   p = Prob.PENOPT.pen; 
   if isempty(x_0)
      x_0 = p.x0(:);
   else
      p.x0 = x_0;
   end
else
   p = [];
end

% Define lower and upper bound arrays for PENSDP
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints

BIG=1E20;
[bl, bu, n, m] = defblbu(Prob, BIG);

% Initial checks on the inputs
%n = Prob.N;

if isempty(x_0)
   x_0 = zeros(n,1);
end

% Safe-guard for x_0
x_0  = max(bl(1:n),min(x_0,bu(1:n)));

% Options
ioptions = DefPar(Prob.PENOPT,'ioptions',-999*ones(12,1));
foptions = DefPar(Prob.PENOPT,'foptions',-999*ones(12,1));

if (length(ioptions) < 12),ioptions=[ ioptions(:) ; -999*ones(12-length(ioptions),1)];end
if (length(foptions) < 12),foptions=[ foptions(:) ; -999*ones(12-length(foptions),1)];end

if ioptions(1) == 0
   ioptions = [];
   foptions = []; % Defaults filled in by MEX;
   
else
   ioptions(1) = 1; % User defines values

   % Check each element and set values
   
   % PBM_MAX_ITER
   if ioptions(2) < 0, ioptions(2) = DefPar(Prob.optParam,'MaxIter',50);end
   % UM_MAX_ITER
   if ioptions(3) < 0, ioptions(3) = DefPar(Prob.optParam,'MinorIter',100);end
   % OUTPUT
   if ioptions(4) < 0, ioptions(4) = 1; end % PRINTLEVEL
   
   if ioptions(5) < 0, ioptions(5) = 0; end % DENSE
   if ioptions(6) < 0, ioptions(6) = 0; end % LS
   if ioptions(7) < 0, ioptions(7) = 0; end % XOUT
   if ioptions(8) < 0, ioptions(8) = 0; end % UOUT
   if ioptions(9) < 0, ioptions(9) = 0; end % NWT_SYS_MODE
   if ioptions(10) < 0, ioptions(10) = 0; end % PREC_TYPE
   if ioptions(11) < 0, ioptions(11) = 0; end % DIMACS
   if ioptions(12) < 0, ioptions(12) = 0; end % TR_MODE
   
   % foptions - floating point values
   
   if foptions(1) < 0, foptions(1) = 1.0; end % U0
   if foptions(2) < 0, foptions(2) = 0.7; end % MU
   if foptions(3) < 0, foptions(3) = 0.1; end % MU2

   %PBM_EPS
   if foptions(4) < 0, foptions(4) = DefPar(Prob.optParam,'eps_f',1e-7);end 
   
   if foptions(5) < 0, foptions(5) = 1e-7;  end % P_EPS
   if foptions(6) < 0, foptions(6) = 1e-14; end % UMIN
   if foptions(7) < 0, foptions(7) = 1e-2;  end % ALPHA
   if foptions(8) < 0, foptions(8) = 1.1;   end % P0

   if foptions(9) < 0, foptions(9) = 0.0;      end % PEN_UP
   if foptions(10) < 0, foptions(10) = 1.0;    end % ALPHA_UP
   if foptions(11) < 0, foptions(11) = 1.0e-7; end % PRECISION_2
   if foptions(12) < 0, foptions(12) = 5.0e-2; end % CG_TOL_DIR

end

% Number of iterations may NOT be less than 2. It causes strange errors in
% the PENSDP solver
if(length(ioptions) >= 2)
  ioptions(2) = max(ioptions(2), 2);
end

%ioptions,foptions
%pause

[mA,nA] = size(Prob.A);
% Check linear part
c = Prob.QP.c(:);
if isempty(c) 
   c = zeros(n,1);
end

if isempty(p)
   
   if ~isempty(Prob.A)
      if nA~=n, error('Linear constraints A MUST have n columns!'); end 
      if mA~=m, error('Linear constraints A MUST have m rows!'); end 
   end
else
end

% Determine type of problem
if isempty(c) | all(c==0)
   Result.f_0=0;
else
   Result.f_0=c(1:n)'*x_0;
end

if isempty(p)
   ixU = find(~isinf(Prob.x_U) & Prob.x_U <  BIG); 
   ixL = find(~isinf(Prob.x_L) & Prob.x_L > -BIG); 
   
   ibU = find(~isinf(Prob.b_U) & Prob.b_U <  BIG); 
   ibL = find(~isinf(Prob.b_L) & Prob.b_L > -BIG); 
   
   % Total number of linear constraints
   % p.ci
   Arhs = full([Prob.b_U(ibU);-Prob.b_L(ibL);Prob.x_U(ixU);-Prob.x_L(ixL)]);
   
   diag = speye(n);
   A = [Prob.A(ibU,:);-Prob.A(ibL,:);diag(ixU,:);-diag(ixL,:)];
   
   SDP = DefPar(Prob.PENOPT,'SDP',[]);

   [x_k, info, Obj, v_k, iresults, fresults] = ...
       pensdp(c, A, Arhs, x_0, SDP, ioptions, foptions, PriLev);
   
else
   p.x0 = x_0;
   [x_k, Obj, v_k, iresults, fresults, info] = pen(p);
end

Inform = info;
Iter = iresults(1);

switch(Inform)
   case 4    , ExitFlag=1;    % Too many iterations
   case {6,7}, ExitFlag = -1; % Fatal errors
   case 5    , ExitFlag = 10; % Input error 
   otherwise , ExitFlag=0;
end

switch(Inform)
   case 0 , Text = 'Solution obtained.';
   case 1 , Text = 'Cholesky factorization of Hessian failed. The result may still be useful.';
   case 2 , Text = 'No progress in objective value, problem probably infeasible.';
   case 3 , Text = 'Linesearch failed. The result may still be useful.';
   case 4 , Text = 'Maximum iteration limit exceeded. The result may still be useful.';
   case 5 , Text = 'Wrong input parameters (ioptions,foptions).';
   case 6 , Text = 'Memory error.';
   case 7 , Text = 'Unknown error, please contact support@tomopt.com';
    otherwise, Text = 'NOTE: UNKNOWN PENBMI Inform value.';
end

Result.ExitText = Text;

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = v_k;

Result.g_k=c;

Result.FuncEv    = Iter;
Result.GradEv    = Iter;
Result.ConstrEv  = Iter;
Result.Iter      = iresults(1);
Result.MinorIter = iresults(2);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

Result.CPUtime   = iresults(4);
Result.REALtime  = iresults(4);

Result.PENOPT.iresults = iresults;
Result.PENOPT.fresults = fresults;
Result.PENOPT.info     = info;


if mA > 0
   Ax = Prob.A * x_k;
else
   Ax = [];
end

% Compute Result.xState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [],Prob.optParam.xTol, ...
   Prob.optParam.bTol, [], bl, bu); 

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nPENSDP solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('PENSDP: Inform = %2d, ',Inform)
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
% 020703 hkh  Written
% 020704 hkh  Generalized to handle PENSDP.pen format and Tomlab QP format
% 020722 ango Added handling of LMI(:,:).Q format of LMI constraints
% 020722 ango Modified for new version of pensdp_w32.lib
% 020725 ango More adaptation towards new version
% 020826 hkh  Revision
% 020826 ango Safe-guarded against empty c-vector
% 020828 ango Help section revised
% 020902 ango Implement ioptions,foptions handling
% 030127 ango Changed LMI(i).Q0 --> LMI(i,1).Q0
% 030711 ango Fixed LMI handling to conform with PENBMI, more efficient and stable.
% 031107 ango PENSDP 1.1 compatible, length of foptions, error codes.
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040803 med  Added pragmas for MATLAB Compiler
% 041203 hkh  Revise calls to defblbu and StateDef
% 041222 med  Safeguard added for x_0
% 041214 frhe Changed to handle the new LMI format
% 041223 frhe The options now work as intended. MaxIter safed to at least 2
% 060818 med  isnan checks removed for x_0
% 070628 ango 2.2 updates to help text
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
