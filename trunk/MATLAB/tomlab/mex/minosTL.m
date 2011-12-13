% TOMLAB MINOS NLP Solver
%
% function Result = minosTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format.
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% c_L, c_U  Bounds on nonlinear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% PriLevOpt Print level.
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
%           1E6 is 10MB of memory. Default n+m1+m2 (number of variables + 
%           linear + nonlinear constraints).
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
% Define: nnL = max(nnObj,nnJac))
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        0         11111   JFLXB:
%                                            Jac, fCon, lambda, x, B=LU stats
% --- Frequencies I
% 5.  PRINT FREQUENCY               0        100
% 6.  SUMMARY FREQUENCY             0        100
% 7.  SOLUTION YES/NO               0        1         1       1 = YES; 0 = NO
% 8.  SUPPRESS PARAMETERS           0        0         1       1 = True
%
% --- Convergence Tolerances
% 9.  ROW TOLERANCE                 >0       1E-6
% 10. OPTIMALITY TOLERANCE          >0       max(1E-6,(10eps_R)^0.5)=1.73E-6
% 11. FEASIBILITY TOLERANCE         >0       1E-6
%
% --- Derivative checking
% 13. VERIFY LEVEL                  -1       -1        3       {-1,0,1,2,3}
% 14. START OBJECTIVE CHECK AT COL  0        1         nnObj
% 15. STOP OBJECTIVE CHECK AT COL   0        nnObj     nnObj
% 16. START CONSTRAINT CHECK AT COL 0        1         nnJac
% 17. STOP CONSTRAINT CHECK AT COL  0        nnJac     nnJac
%
% --- Scaling
% 18. SCALE OPTION                  0        1 or 2    2       2 if LP,1 if NLP
%     See the TOMLAB /MINOS manual for more information.
% 19. SCALE TOLERANCE               >0       0.9       <1
% 20. SCALE PRINT                   0        0         1       1 = True
% 21. CRASH TOLERANCE               0        0.1       <1
%
% --- Other Tolerances
% 22. LINESEARCH TOLERANCE          >0       0.1       <1
%
% --- LU I
% 23. LU FACTORIZATION TOLERANCE    1        100 or 5         100 if LP
% 24. LU UPDATE TOLERANCE           1        10  or 5          10  if LP
% 25  LU SWAP TOLERANCE             >0       1.22E-4           eps^(1/4)
% 26. LU SINGULARITY TOLERANCE      >0       3.25E-11          eps^(0.67)
%
% --- LP or LC subproblems
% 27. PIVOT TOLERANCE               >0       3.25E-11          eps^(0.67)
% 28. CRASH OPTION                  0        3         3       {0,1,2,3}
% 29. WEIGHT ON LINEAR OBJECTIVE    0.0      0.0               during Phase 1
% 30. ITERATION LIMIT               0        3(m+m3) + 10nnL
%     m3=1 if length(Prob.QP.c) > 0, otherwise m3=0
%     TOMLAB default: max(10000,3(m+m3) + 10nnL)
% 31. PARTIAL PRICE                 1        10 or 1           10 if LP
%
% --- SLC method
% 32. MAXIMIZE                      0        0         1       1=maximize
% 33. LAGRANGIAN                    0        1         1       1=YES, 0=NO
% 34. PENALTY PARAMETER             0.0      1.0
% 35. MAJOR ITERATIONS LIMIT        >0       50
% 36. MINOR ITERATIONS LIMIT        >0       40
% 37. MAJOR DAMPING PARAMETER       >0       2.0
% 38. MINOR DAMPING PARAMETER       >0       2.0
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     Is always set by minosTL dependent on Prob.ConsDiff, Prob.NumDiff
% 40. RADIUS OF CONVERGENCE         0.0      0.01
% 41. FUNCTION PRECISION            >0       3.00E-13          eps^0.8=eps_R
% 42. DIFFERENCE INTERVAL           >0       5.48E-8           eps^0.4
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.69E-5           eps^{0.8/3}
% 44. COMPLETION                    0      1 LC, 0 NC  1     0=PARTIAL 1=FULL
% 45. UNBOUNDED STEP SIZE           >0       1E10
% 46. UNBOUNDED OBJECTIVE           >0       1E20
%
% --- Hessian approximation
% 47. HESSIAN DIMENSION             1        50        1+nnL
% 48. SUPERBASICS LIMIT             1        50        1+nnL
%     TOMLAB default (to avoid termination with Superbasics Limit too small):
%     If n <= 5000: max(50,n+1)
%     If n >  5000: max(500,n+200-size(A,1)-length(c_L))
%     Avoid setting REDUCED HESSIAN (number of columns in reduced Hessian).
%     It will then be set to the same value as the SUPERBASICS LIMIT by MINOS
%
% --- Frequencies II
% 51. CHECK FREQUENCY               >0       60
% 52. EXPAND FREQUENCY              >0       10000
% 53. FACTORIZATION FREQUENCY       >0       50
%
% --- LU II
% 63. LU PARTIAL  PIVOTING          0        0         3  0=partial
%     or LU COMPLETE PIVOTING                             1=complete
%     or LU ROOK PIVOTING                                 2=rook
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
% x_0      Initial  solution vector.
% g_k      Gradient of the function.
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from minos.m (similar to TOMLAB).
% Inform   MINOS information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
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
% For a problem description, see conAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written July 10, 2000. Last modified June 4, 2008.

% -----------------------------------------------
% Experimental parameters, only available using SPECS file:
% -----------------------------------------------
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
% xx. MULTIPLE PRICE                1        1                 >1 often bad
% xx. LU DENSITY TOLERANCE          >0       0.5         1

function Result = minosTL(Prob)

%#function nlp_cdcS nlp_fg lp_f lp_g lp_H

if nargin < 1, error('minosTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 3; % NLP (CON) solver

if isempty(Prob.NumDiff)
   ObjDers = 1;
elseif Prob.NumDiff==6
   ObjDers = 0;
else
   ObjDers = 1;
end
if Prob.mNonLin > 0
   if isempty(Prob.ConsDiff)
      ConsDers = 1;
   elseif Prob.ConsDiff==6
      ConsDers = 0;
   else
      ConsDers = 1;
   end
else
   ConsDers = 0;
end

Prob = iniSolve(Prob,3,ObjDers,ConsDers);

Result=ResultDef(Prob);
Result.Solver='MINOS';
Result.SolverAlgorithm='MINOS 5.51 NLP code';

% Define lower and upper bound arrays for MINOS
%
% Inf are changed to BIG (=1E10), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints

BIG=1E10;

% bl/bu has order (x, nonlinear constraints, linear constraints)
[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 0);

PriLev=Prob.PriLevOpt;

% Initial checks on the inputs
ObjAdd  = 0;
m       = m1+m2;

if isfield(Prob.SOL,'moremem')
   moremem = Prob.SOL.moremem;
else
   moremem = [];
end
if isempty(moremem), moremem = n+m; end

% Check if any linear part
if isfield(Prob.QP,'c')
   c  = Prob.QP.c;
   m3 = length(c) > 0;
else
   c  = [];
   m3 = 0;
end

% Check on Warm start, then set hs,xs,nS

if Prob.WarmStart
   % Warm start for SOL solver
   Warm = 1;
   xs   = Prob.SOL.xs;
   hs   = Prob.SOL.hs;
   nS   = Prob.SOL.nS;
   x_0  = xs(1:n);
else
   % Initial values
   x_0  = Prob.x_0(:);
   Warm = 0;
   nS   = Prob.SOL.nS;
   hs   = Prob.SOL.hs;
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
   xs   = x_0;
end
% Avoid putting c both in objective, and calling lp_f to get c'x
if strcmpi(Prob.FUNCS.f,'lp_f')
   nnObj = 0;
else
   nnObj = n; % number of nonlinear variables
end

if m2 > 0
   % Determine the sparse problem structure

   % number of nonlinear constraints   (nnCon)
   nnCon=m2;
   if isempty(Prob.ConsPattern)
      gJac=ones(m2,n);
      nnJac=n;
   else
      gJac=Prob.ConsPattern;
      [ix,iy]=find(gJac);
      % [i,j]=find(Prob.ConsPattern);

      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(gJac),ix,iy);

      % Number of nonlinear jacobian vars (nnJac)
      nnJac=size(gJac,2);
   end
   % Must call function, in case global values are needed for constraints
   % Because minos calls constraints first, then functions
   Result.f_0 = nlp_f(x_0, Prob);
   % Sometimes maybe also a call to the gradient might be needed.
   % Do not make such a call now, it is highly unlikely.
else
   nnCon      = 0;
   nnJac      = 0;
   gJac       = [];
   if ~isempty(x_0) & ~isempty(c)
      Result.f_0 = c'*x_0;
   else
      %Result.f_0 = 0;
      Result.f_0 = nlp_f(x_0, Prob);
   end
end

% Must change sign, and interchange the constraint bounds, 
% i.e. set correct bounds on slack variables

v=-bl(n+1:n+m);

bl(n+1:n+m)=-bu(n+1:n+m);
bu(n+1:n+m)=v;

% NOT MORE!!! Also move nonlinear constraints BEFORE linear constraints
%bl(n+1:n+m)=-bu([n+m1+1:n+m,n+1:n+m1]);
%bu(n+1:n+m)=v([m1+1:m,1:m1]);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end 
end 

% Set up the constraint matrix A 
if isempty(gJac)
   A=sparse(Prob.A);
else
   if isempty(Prob.A)
      A = sparse(gJac);
   else
      if nnJac < n
         A = sparse([[gJac,zeros(m2,n-nnJac)];Prob.A]);
      else
         A = sparse([gJac;Prob.A]);
      end
   end
end

if m==0
   % Construct one linear constraint from some bounded variable
   i = find(bu < BIG);
   if isempty(i)
      i = find(bl > -BIG);
   end
   if isempty(i), i=1; end
   i = i(1);

   bl(n+1,1)=-min(bu(i),BIG);
   bu(n+1,1)=-max(bl(i),-BIG);
   m=1;
   A = sparse(zeros(1,n));
   A(1,i)=1;
end

if m3 > 0
   A = [A;c(:)'];
   iObj=m+1;  % The constraint row with linear objective term c is put last
   % Add the bounds for the objective in bl and bu
   bl=[bl;-BIG];
   bu=[bu; BIG];
   % Must set Prob.QP.c as empty, if the problem is a pure QP or LP.
   Prob.QP.c=[];
else
   iObj=0;
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

% Lagrange multipliers, known info
pi(1:m+m3)=0;

% No QP (No check if it happens to be a QP or LP), H set as empty
H=[];

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('minos',3,...
%         nnObj, nnJac, size(A,1), Prob);

optPar    = Prob.SOL.optPar(:)';
optPar(length(optPar)+1:Prob.SOL.optParN)=-999;

if optPar(30) <= 0
   % Tomlab default: max(2000,3(m+m3) + 10nnL)
   nnL = max(nnObj,nnJac);
   optPar(30) = max(2000,3*(m+m3)+10*nnL);
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
if optPar(48) < 0 & (n-m1-m2) > 45
   % Increase number of superbasics
   if n > 5000
      optPar(48) = max(500,n-m1-m2+200);
   else
      optPar(48) = max(50,n+1);
   end
end
% Avoid derivative check
if optPar(13) < 0
   optPar(13)=-1;
end

[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     minos( H, A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
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
if Inform == 1
   %z = xs(1:n);
   % Adjust infeasible solution inside lower and upper variable bound
   xs(1:n) = max(bl(1:n),min(bu(1:n),xs(1:n)));
   %ixx = find(xs(1:n) ~= z);
   %fprintf('Adjust %d values\n',length(ixx))
end

Result.f_k=Obj;
Result.g_k=gObj;
Result.c_k=fCon;
Result.x_k=xs(1:n);
Result.x_0=x_0;
if m1 + m2 == 0
   Result.v_k=rc(1:n);
else
   Result.v_k=[rc(1:n);pi(1:m)];
end


% Saved for a warm start
% The dummy objective row is last in xs, value -Obj, last in hs is 3 (basic)
Result.SOL.xs=xs;
Result.SOL.hs=hs;
Result.SOL.nS=nS;
Result.SOL.nInf=nInf;
Result.SOL.sInf=sInf;
Result.SOL.optPar=optPar;

if ~isempty(gCon)
   [ix,iy] = find(gJac);
   % TOMLAB now has the format one constraint / row, same as SOL solvers
   Result.cJac=sparse(ix,iy,gCon,size(gJac,1),size(gJac,2));
end

Result.FuncEv    = [];% iwCount(3);
Result.GradEv    = [];% sum(iwCount(4:6));
Result.ConstrEv  = [];% iwCount(7);
Result.Iter      = iwCount(1);
Result.MinorIter = iwCount(2);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

optParam = Prob.optParam;
if m1 > 0
  Result = StateDef(Result, xs(1:n), A(m2+1:m2+m1,:)*xs(1:n), fCon, ...
                 optParam.xTol, optParam.bTol, optParam.cTol, ...
                 [bl(1:n);-bu(n+1:n+m1+m2)], [bu(1:n);-bl(n+1:n+m1+m2)],0); 
else
  Result = StateDef(Result, xs(1:n), [], fCon, ...
                 optParam.xTol, optParam.bTol, optParam.cTol, ...
                 [bl(1:n);-bu(n+1:n+m1+m2)], [bu(1:n);-bl(n+1:n+m1+m2)],0); 
end

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
   fprintf('\nMINOS solving Problem %d:\n',Prob.P);
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
   fprintf('fObj and gObj evaluations%7d %d %d %d\n',iwCount(3:6));
   if m2 > 0
      fprintf('fCon and gCon evaluations%7d %d %d %d\n',iwCount(7:10));
   end
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
elseif Prob.optParam.IterPrint 
   fprintf('MINOS: Inform =%2d, ',Inform)
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000710 hkh  New minosMex interface.
% 000828 hkh  If file name given, and print level in optPar(1) is 0, set to max
% 000916 hkh  Return string ExitText with interpretation of Inform flag
% 001005 hkh  Must set Prob.QP.c as empty, if called with a pure LP or QP
% 010901 hkh  Give minos correct number of nonlinear Jacobian variables
%             using size of Prob.ConsPattern
% 010903 hkh  Add complete pivoting option. Allow cold start change of hs,ns
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 020210 hkh  Send linear index for sparse gJac subscripts in Prob.ConsIdx
% 020304 hkh  Set optPar(39) DERLVL dependent on Prob.ConsDiff, Prob.NumDiff
% 020417 hkh  Send optPar back in Result.SOL
% 030406 hkh  Adjust infeasible solution inside lower and upper variable bound
% 030406 hkh  Print Inform if Prob.optParam.IterPrint = 1, and PriLev == 0
% 031118 hkh  Must change BIG from 1E20 to 1E10, otherwise incorrect solution
% 031120 ango Help updated for optPar(63), for MINOS 5.51
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040602 med  Help fix PriLev to PriLevOpt
% 040728 med  Pragmas added for MATLAB Compiler
% 040929 hkh  Default line search tolerance is 0.1, not 0.9, for minos
% 041202 hkh  If internal differentiation, different call to iniSolve
% 041202 hkh  Revise call to defblbu and StateDef, use different Order=0
% 041222 med  Safeguard added for x_0
% 050221 med  Help about scaling updated
% 050606 hkh  Add new parameters AIJ CONVERGENCE, SUBSPACE, LU WRAP
% 050606 hkh  Change defaults for #6,42,43
% 050606 hkh  Increase SUPERBASICS, optPar(48) to avoid too low value
% 050614 med  Updated help for optPar vector
% 050614 med  Removed optPar 3 and 4, 54-62
% 050614 hkh  Removed SOLSet call, set and comment optPar 30,39,48
% 050614 hkh  Use hs to compute StateDef fields
% 050617 med  optPar defaults 45 and 46 swapped
% 050725 med  bl, bu safe-guarded to be a column vector
% 050829 med  FuncEv, GradEv, ConstrEv set to []
% 050829 med  Added correct StateDef
% 051216 med  optPar(13) default to -1 (no derivative check)
% 051216 med  Inform 6, 7 and 8 text updated
% 060818 med  isnan checks removed for x_0
% 061025 med  Linear part of objective fixed for nonlinear problems
% 061224 med  moremem parameter added
% 080604 hkh  Safeguard call to iniSolve for empty NumDiff and ConsDiff
