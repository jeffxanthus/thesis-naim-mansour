% TOMLAB MINOS QP Solver
%
% function Result = minosqpTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj. nnObj < n is OK.
% PriLevOpt Print Level.
% WarmStart If true, use warm start, otherwise cold start. Use WarmDefSOL
%           to set the proper parameters.
%
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution and slacks from previous run.
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%                 2=superbasic (between bounds), 3=basic (between bounds)
% nS        Number of superbasics from previous run.
% moremem   Add more memory if MINOS stops with not enough storage message.
%           1E6 is 10MB of memory. Default n+m (number of variables +
%           linear constraints).
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined
%           by SPECS parameters or optPar parameters.
% SummFile  Name of SOL Summary File.
% optPar    Elements > -999 take effect.
%
% -----------------------------------------------
% How optPar is used for setting SPECS parameters:
% -----------------------------------------------
%
% optPar  Structure with optimization parameters.
%
% MINOS keywords in optPar(#):
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        0         1    0=brief 1=LU stats
%
% --- Frequencies I
% 5.  PRINT FREQUENCY               0        100
% 6.  SUMMARY FREQUENCY             0        1
%
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS PARAMETERS           0        0         1       1 = True
%
% --- Convergence Tolerances
% 10. OPTIMALITY TOLERANCE          >0       max(1E-6,(10eps_R)^0.5) = 1.73E-6
% 11. FEASIBILITY TOLERANCE         >0       1E-6
%
% --- Derivative checking
% 13. VERIFY LEVEL                  -1       -1        3       {-1,0,1,2,3}
% 14. START OBJECTIVE CHECK AT COL  0        1         nnObj
% 15. STOP OBJECTIVE CHECK AT COL   0        nnObj     nnObj
%
% --- Scaling
% 18. SCALE OPTION                  0        1         2
%     See the TOMLAB /MINOS manual for more information.
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
% 21. CRASH TOLERANCE               0        0.1       <1
%
% --- Other Tolerances
% 22. LINESEARCH TOLERANCE          >0       0.1       <1
%
% --- LU I
% 23. LU FACTOR TOLERANCE           1        5.0
% 24. LU UPDATE TOLERANCE           1        5.0
% 25  LU SWAP TOLERANCE             >0       1.22E-4           eps^(1/4)
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
%
% --- LP parameters
% 27. PIVOT TOLERANCE               >0       3.25E-11          eps^(0.67)
% 28. CRASH OPTION                  0        3         3       {0,1,2,3}
% 29. WEIGHT ON LINEAR OBJECTIVE    0.0      0.0               during Phase 1
% 30. ITERATION LIMIT               0        3(m+m3) + 10nnL
%     m3=1 if length(Prob.QP.c) > 0, otherwise m3=0
%     Tomlab default: max(10000,3(m+m3) + 10nnObj)
% 31. PARTIAL PRICE                 1        1
% 32. MAXIMIZE                      0        0         1       1=maximize
%
% --- Reduced-gradient method
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     Is always set by minosqpTL to 3
% 41. FUNCTION PRECISION            >0       3.0E-13           eps^0.8=eps_R
% 42. DIFFERENCE INTERVAL           >0       5.48E-8           eps^0.4
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.70E-5           eps^{0.8/3}
% 44. COMPLETION                    0      1 LC, 0 NC  1     0=PARTIAL 1=FULL
% 45. UNBOUNDED STEP SIZE           >0       1E10
% 46. UNBOUNDED OBJECTIVE           >0       1E20
%
% --- Hessian approximation
% 47. HESSIAN DIMENSION             1        50        1+nnObj
% 48. SUPERBASICS LIMIT             1        50        1+nnObj
%     TOMLAB default (to avoid termination with Superbasics Limit too small):
%     If n <= 5000: max(50,n+1)
%     If n >  5000: max(500,n+200-size(A,1)-length(c_L))
%     Avoid setting REDUCED HESSIAN (number of columns in reduced Hessian).
%     It will then be set to the same value as the SUPERBASICS LIMIT by MIN
%
% --- Frequencies II
% 51. CHECK FREQUENCY               >0       60
% 52. EXPAND FREQUENCY              >0       10000
% 53. FACTORIZATION FREQUENCY       >0       50
%
% --- LU II
% 63. LU COMPLETE PIVOTING          0        0         2    2=rook, 1=complete,0=partial
%     or LU PARTIAL PIVOTING
%     or LU ROOK PIVOTING
%
% --- Additional parameters
% 67. AIJ TOLERANCE                 0     1E-10
%     Elements |a(i,j)| < AIJ TOLERANCE are set as 0
% 70  SUBSPACE                     >0      0.5         1  Subspace tolerance
%     Convergence tolerance in current subspace before consider moving off
%     another constraint
%
% -----------------------------------------------------------------------
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
% ExitFlag Exit status from minos.m (similar to TOMLAB).
% Inform   MINOS information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (minos).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% xs       Solution and slack variables.
% hs       State for variables and slacks in xs.
% nS       # of superbasics.
% nInf     # of infeasibilities.
% sInf     Sum of infeasibilities.
%
% -----------------------------------------------------------------------
%
% For a problem description, see qpAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written July 10, 2000.  Last modified Jul 17, 2009.

function Result = minosqpTL(Prob)

if nargin < 1, error('minosqpTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 3;

Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='MINOS';
Result.SolverAlgorithm='MINOS 5.51 QP (NLP) code';

PriLev=Prob.PriLevOpt;

BIG=1E10;

% m2 = 0 for QPs
[bl, bu, n, m, m2] = defblbu(Prob, BIG, 1);
A = Prob.A;
ObjAdd  = 0;

% Initial checks on the inputs

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isempty(moremem), moremem = n+m; end

% Check if any linear part
c  = Prob.QP.c;
m3 = ~isempty(c);
nb = n+m+m3;

% Check on Warm start, then set hs,xs,nS

if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs;
   hs   = Prob.SOL.hs;
   nS   = Prob.SOL.nS;
   x_0  = xs(1:n);
else
   Warm = 0;
   nS   = Prob.SOL.nS;
   hs   = Prob.SOL.hs;
   % Initial values
   x_0  = Prob.x_0(:);
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
   xs   = x_0;
end

nnObj  = size(Prob.QP.F,1); % number of nonlinear variables

if isempty(A)
   % Construct one linear constraint from some bounded variable when m==0
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);
   bl=[bl;-min(bu(i),BIG)];
   bu=[bu;-max(bl(i),-BIG)];
   A = sparse(zeros(1,n));
   A(1,i)=1;
   if length(hs) < n+1
      hs=[hs;0];
   end
end

if m3 > 0
   A = [sparse(A);c'];
   iObj=m+1;  % The constraint row with linear objective term c is put last
   % Add the bounds for the objective in bl and bu
   bl=[bl;-BIG];
   bu=[bu; BIG];
else
   A=sparse(A);
   iObj=0;
end

if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Non linear variables %d. ',nnObj)
      fprintf('Total number of constraints %d.\n',m)
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
      mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 10, pause; end
   end
end


% Must change sign, and interchange the constraint bounds,
% i.e. set correct bounds on slack variables

v=-bl(n+1:n+m);
bl(n+1:n+m)=-bu(n+1:n+m);
bu(n+1:n+m)=v;

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

if isempty(c)
   if nnObj==0
      Result.f_0=0;
   else
      Result.f_0=0.5*(x_0(1:nnObj)'*Prob.QP.F*x_0(1:nnObj));
   end
elseif nnObj==0
   Result.f_0=c(1:n)'*x_0(1:n);
else
   Result.f_0=0.5*(x_0(1:nnObj)'*Prob.QP.F*x_0(1:nnObj)) + c(1:n)'*x_0(1:n);
end

pi    = [];
nnObj = size(Prob.QP.F,1);

optPar    = Prob.SOL.optPar(:)';
optPar(length(optPar)+1:Prob.SOL.optParN)=-999;

if optPar(30) <= 0
   % Tomlab default: max(2000,3(m+m3) + 10nnObj)
   optPar(30) = max(2000,3*(m+m3)+10*nnObj);
end

if Prob.NumDiff == 6
   if Prob.ConsDiff == 6
      optPar(39)=0;
   else
      optPar(39)=2;
   end
elseif Prob.ConsDiff == 6
   optPar(39)=1;
else
   optPar(39)=3;
end
% Avoid termination due to 'The Superbasics limit is too small'
if optPar(48) < 0 & (n-m-m2) > 45
   % Increase number of superbasics
   if n > 5000
      optPar(48) = max(500,n-m-m2+200);
   else
      optPar(48) = max(50,n+1);
   end
end
% Avoid derivative check
if optPar(13) < 0
   optPar(13)=-1;
end

[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj] = ...
     minos( Prob.QP.F, A, bl, bu, 0, nnObj, 0, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            Prob.SOL.SpecsFile, Prob.SOL.PrintFile, Prob.SOL.SummFile, ...
            PriLev, ObjAdd, moremem, Prob.Name );

switch Inform
   case 3
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 1
     ExitFlag=4;  % Infeasible
   case {10,22}
     ExitFlag=3;  % Rank problem
   case {5,7,8,20,21,30,32,42,43,44}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k   = Obj;
Result.x_k   = xs(1:n);
Result.x_0   = x_0;
Result.v_k   = [rc(1:n);pi(1:m)];

if ~isempty(c)
   if nnObj > 0
      Result.g_k=Prob.QP.F*xs(1:n)+c;
   else
      Result.g_k=c;
   end
elseif nnObj > 0
   Result.g_k=Prob.QP.F*xs(1:n);
else
   Result.g_k=zeros(n,1);
end

% Saved for warm start
% The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
Result.SOL.xs=xs;
Result.SOL.hs=hs;
Result.SOL.nS=nS;
Result.SOL.nInf=nInf;
Result.SOL.sInf=sInf;

Iter=iwCount(2);

Result.FuncEv    = Iter;
Result.GradEv    = 0;
Result.ConstrEv  = Iter;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

Result = StateDef(Result, xs(1:n), A(1:m,:)*xs(1:n), [], ...
                 Prob.optParam.xTol, Prob.optParam.bTol, [], ...
                 [bl(1:n);-bu(n+1:n+m)], [bu(1:n);-bl(n+1:n+m)],1);

switch Inform
   case 0
      Text = 'Optimal solution found';
   case 1
      Text = 'The problem is infeasible';
   case 2
      Text = 'The problem is unbounded (or badly scaled)';
   case 3
      Text = 'Too many iterations';
   case 4
      Text = str2mat('Apparent stall.  The solution has not changed' ...
             ,'for a large number of iterations (e.g. 1000).');
   case 5
      Text = 'The Superbasics limit is too small';
   case 6
      Text = 'User requested termination (by returning bad value)';
   case 7
      Text = 'Gradient seems to be giving incorrect derivatives';
   case 8
      Text = 'Jacobian seems to be giving incorrect derivatives';
   case 9
      Text = 'The current point cannot be improved';
   case 10
      Text = str2mat('Numerical error in trying to satisfy the linear ' ...
             ,'constraints (or the linearized nonlinear constraints)' ...
             ,'The basis is very ill-conditioned.');
   case 11
      Text = 'Cannot find a superbasic to replace a basic variable';
   case 12
      Text = str2mat('Basis factorization requested twice in a row' ...
             ,'Like case Inform = 9. Possibly convergence?');
   case 13
      Text = str2mat('Near-optimal solution found' ...
             ,'Should probably be treated as inform = 9');
   case 20
      Text = 'Not enough storage for the basis factorization';
   case 21
      Text = 'Error in basis package';
   case 22
      Text = str2mat('The basis is singular after several attempts' ...
            ,'to factorize it (and add slacks where necessary)');
   case 30
      Text = str2mat('An OLD BASIS file had dimensions that did' ...
             ,'not match the current problem');
   case 32
      Text = 'System error. Wrong number of basic variables';
   case 40
      Text = 'Fatal errors in the MPS file';
   case 41
      Text = 'Not enough storage to read the MPS file';
   case 42
      Text = 'Not enough storage to solve the problem';
   otherwise
      Text = 'NOTE: UNKNOWN MINOS Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nMINOS QP solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('MINOS: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Major iterations%7d. ',iwCount(1));
   fprintf('Minor iterations%7d. ',iwCount(2));
   fprintf('\n');
   fprintf('nInf            %7d. ',nInf);
   fprintf('nS              %7d. ',nS);
   fprintf('sInf %14.7e\n',sInf);

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(xs(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector hs for x and constraints = \n');
      xprinti(hs(1:min(length(hs),MAX_x)),'hs: ');
   end
   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (pi) = \n');
      xprinte(pi(1:min(length(pi),MAX_c)),'pi:');
      fprintf('Reduced costs rc: Last %d elements should be -v_k\n',m);
      xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
   end
end

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 000714 hkh  New minosMex QP interface.
% 000828 hkh  If file name given, and print level in optPar(1) is 0, set to max
% 000916 hkh  Return string ExitText with interpretation of Inform flag
% 001005 hkh  Revision for new minos.dll
% 010903 hkh  Add complete pivoting option. Allow cold start change of hs,ns
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 031118 hkh  Must change BIG from 1E20 to 1E10, otherwise incorrect solution
% 031120 ango Help for optPar(63) changed, MINOS 5.51
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040413 hkh  Make full call to iniSolve
% 040602 med  Help fix PriLev to PriLevOpt
% 041004 med  Linesearch tolerance 0.1 default
% 041202 hkh  Revise call to defblbu and StateDef
% 041222 med  Safeguard added for x_0
% 050221 med  Help about scaling updated
% 050614 med  Updated help for optPar vector
% 050614 med  Removed optPar 3 and 4, and 54-62
% 050614 hkh  Removed SOLSet call, set and comment optPar 30,39,48
% 050614 hkh  Use hs to compute StateDef fields
% 050617 med  optPar defaults 45 and 46 swapped
% 050725 med  Modified Prob.A no longer passed to Result
% 051216 med  optPar(13) default to -1 (no derivative check)
% 051216 med  Inform 6, 7 and 8 text updated
% 060818 med  isnan checks removed for x_0
% 061224 med  moremem parameter added
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
% 090717 med  f_0 calculation updated
