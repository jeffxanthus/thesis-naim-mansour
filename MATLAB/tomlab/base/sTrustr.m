% sTrustr.m:  Structured Trust Region algorithm
% Ref. A. R. Conn, Nick Gould, A. Sartenaer, Ph. L. Toint,
% "Convergence Properties of Minimization Algorithms For Convex
% Constraints Using A Structured Trust Region", July 4, 1995
%
% sTrustr implements an algorithm for solving
% linearly constrained  minimization problems
% It handles nonlinear constraints, but might fail if the constraints are
% nonconvex.
% Run e.g. conSolve or nlpSolve for nonlinearly constrained problems.
%
%	min      f(x) = sum (i=1:M) f_i (x), i.e. f is partially separable
%
%       s/t       x_L <=   x  <= x_U
%                 b_L <= A x  <= b_U
%                 c_L <= c(x) <= c_U.
%
% function Result = sTrustr(Prob, varargin)
%
% sTrustr runs different methods to obtain the gradient g and the Hessian H,
% dependent on the parameters Prob.Solver.Alg, Prob.NumDiff, Prob.ADObj
% and if a user supplied routine for the Hessian, stored in
% Prob.FUNCS.H is available.
% Solver.Alg=1 gives quasi-Newton BFGS and Solver.Alg=2 gives DFP
%
% The table gives the different possibilities
%
% Solver.Alg  NumDiff    ADObj isempty(FUNCS.H)  Hessian computation
%     0          0         0           0         Analytic Hessian
%     0          0         0           1         Numerical differences H
%     0         >0         0         any         Numerical differences g,H
%     0         <0         0         any         Numerical differences H
%     0        any        -1         any         Automatic differentiation
%     1          0         0         any         BFGS
%     2          0         0         any         DFP
%     1        ~=0         0         any         BFGS, numerical gradient g
%     2        ~=0         0         any         DFP,  numerical gradient g
%     1        any         1         any         BFGS, automatic diff gradient
%     2        any        -1         any         DFP,  automatic diff gradient
%
% INPUT PARAMETERS
%
%   Use conAssign.m (or probAssign) to initialize the Prob structure in the
%   TOMLAB Quick format.  Fields used in structure Prob:
%
%
%   PartSep Field for partially separable function (psf) information:
%   PartSep.pSepFunc Number of elements M in the psf
%   PartSep.index    Index for the function to compute
%
%   x_0      Starting point
%   x_L      Lower bounds for x
%   x_U      Upper bounds for x
%   b_L:     Lower bounds for linear constraints
%   b_U:     Upper bounds for linear constraints
%   A:       Linear constraint matrix
%   c_L:     Lower bounds for nonlinear constraints
%   c_U:     Upper bounds for nonlinear constraints
%   ConsDiff Differentiation method for the constraint Jacobian
%            0 = analytic, 1-5 different numerical methods
%   SolverQP Name of the solver used for QP subproblems. If empty,
%            picked from a list, best available with a license
%   SolverFP Name of the solver used for FP (feasible point) subproblems.
%            If empty, picked from a list, best available with a license
%   f_Low    A lower bound on the optimal function value.
%            Used in convergence tests, f_k(x_k) <= f_Low. Only a feasible
%            point x_k is accepted
%
%
%   optParam Optimization parameters
%   optParam.QN_InitMatrix  Initial Quasi-Newton matrix, if not empty,
%
% optParam structure in Prob. Fields used:
%    IterPrint     Print short information each iteration
%    eps_f         Relative tolerance in f(x)
%    eps_g         Gradient convergence tolerance
%    eps_x         Convergence tolerance in x
%    eps_absf      Lower bound on function value
%    size_x        Approximate size of optimal variable values
%    size_f        Approximate size of optimal function value
%    MaxIter       Maximal number of iterations
%    wait          Pause after printout if true
%    cTol          Constraint violation convergence tolerance
%    bTol          Linear constraint violation convergence tolerance
%    xTol          Variable violation tolerance
%    PriLev        Print level
%    IterPrint     Print short information each iteration
%    QN_InitMatrix Initial Quasi-Newton matrix, if not empty,
%                  Otherwise use identity matrix
%
% The Prob structure could be created in the TOMLAB format with
% call a call to conAssign,
% or in the Init File format, see Users Guide.
%
% ---------------------------------------- extra parameters:
% VARARGIN: Extra user parameters passed to the computation of f(x) etc
%
% OUTPUT PARAMETERS
%
% Result      Structure with results from optimization
%    x_k      Optimal point
%    v_k      Lagrange multipliers
%    f_k      Function value at optimum
%    g_k      Gradient vector at optimum
%    x_0      Starting value vector
%    B_k      The Quasi Newton matrix if Result.Solver.Alg > 0
%    H_k      The Hessian matrix if Result.Solver.Alg == 0
%    xState   Variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%    Iter     Number of iterations
%    ExitFlag Flag giving exit status
%    ExitText Text string giving ExitFlag and Inform information
%    Inform   Code telling type of convergence
%             1   Iteration points are close.
%             2   Projected gradient small.
%             3   Iteration points are close and projected gradient small.
%             4   Relative function value reduction low for LowIts iterations.
%             5   Iteration points are close and relative function value
%                 reduction low for LowIts iterations.
%             6   Projected gradient small and relative function value
%                 reduction low for LowIts iterations.
%             7   Iteration points are close, projected gradient small and
%                 relative function value reduction low for LowIts iterations.
%             8   Too small trust region
%             9   Trust region small. Iteration points close
%             10  Trust region and projected gradient small.
%             11  Trust region and projected gradient small, iterations close.
%             12  Trust region small, Relative f(x) reduction low.
%             13  Trust region small, Relative f(x) reduction low.
%                 Iteration points are close.
%             14  Trust region small, Relative f(x) reduction low.
%                 Projected gradient small.
%             15  Trust region small, Relative f(x) reduction low.
%                 Iteration points close, Projected gradient small.
%             101 Max no of iterations reached.
%             102 Function value below given estimate.
%             103 Too large trust region, failure
%
% Printing (PriLev = Prob.PriLevOpt):
%
% PriLev: < 0: Totally silent,          0: Error messages,
%         1: Final result,              2: Each iteration, short output
%         3: Each iteration, more info. 4: Info from QP solver. 5: Hessian

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1996-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written Sept 12, 1996.  Last modified Jul 17, 2009.

function Result = sTrustr(Prob, varargin)

%#function lp_f lp_g lp_H

if nargin < 1
   error('sTrustr needs input structure Prob');
end

global MAX_c % Max number of constraints to print

% Calls internal routines itrr, SafeUpdate

solvType=checkType('con');

Prob=ProbCheck(Prob,'sTrustr',solvType);

Prob = iniSolve(Prob,solvType,2,2);

DEBUG=0;     % Debugging flag

% Pick up input variables from Prob structure

% Function names (text strings)
H_k  = [];
B_k  = [];
userH= Prob.FUNCS.H;

[n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

Prob.x_L=x_L;
Prob.x_U=x_U;

lam_L = zeros(n,1);     % Initial zeros for Lagrange multipliers

idx=[1:n];  % No check on which variables to update in BFGS/DFP

% Linear constraints

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

Prob.b_L=b_L;
Prob.b_U=b_U;

% Nonlinear constraints

[m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin{:} );

Prob.c_L=c_L;
Prob.c_U=c_U;

mTot=mA+m;

optParam=Prob.optParam;

eps_g =  optParam.eps_g;    % Gradient convergence tolerance
eps_x =  optParam.eps_x;    % Convergence tolerance in x
MaxIter= optParam.MaxIter;  % Maximal number of iterations
size_x = optParam.size_x;   % Approximate size of optimal variable values
size_f = optParam.size_f;   % Approximate size of optimal function value

fLow   =   Prob.f_Low;      % Lower bound on function value

% If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bounds
xTol   =  optParam.xTol;
bTol   =  optParam.bTol;
cTol   =  optParam.cTol;    % Nonlinear constraint feasibility tolerance
% If (f_km1-f_k) < eps_f * max(f_k,size_f) for LowIts iter, stop
eps_f  = optParam.eps_f;

zTol=[bTol*ones(mA,1);cTol*ones(m,1)];

if MaxIter <= 0, MaxIter=100*n; end

Result=ResultDef(Prob);

PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration

% M is the number of terms in the partially separable functions in f(x)
if isfield(Prob.PartSep,'pSepFunc')
   M = Prob.PartSep.pSepFunc;
else
   M=[];
end
if isempty(M), M=1; end

[alg, SolvAlg, B_k, LowIts]=NameAlg(Prob,M,userH,PriLev);


Result.Solver='sTrustr';
Result.SolverAlgorithm=SolvAlg;

if isempty(Prob.SolverQP)
   Convex = alg == 1 | alg == 2;
   if Convex
      if checkLicense('sqopt')
         SolverQP='sqopt';
      elseif checkLicense('lssol')
         SolverQP='lssol';
      elseif checkLicense('qpopt')
         SolverQP='qpopt';
      elseif checkLicense('bqpd')
         SolverQP='bqpd';
      else
         SolverQP='qld';
      end
   else
      if checkLicense('qpopt')
         SolverQP='qpopt';
      elseif checkLicense('bqpd')
         SolverQP='bqpd';
      else
         SolverQP='qpSolve';
      end
   end
   % SolverQP=GetSolver('qp',Prob.LargeScale,Convex);
   Prob.SolverQP = SolverQP;
else
   SolverQP=Prob.SolverQP;
end

% Setup Structure used in QP calls
ProbQP = CreateProbQP(Prob, 'qp', max(500,3*n), PriLev-3, optParam.wait);

if isempty(Prob.SolverFP)
   if checkLicense('minos')
      SolverFP='lp-minos';
   elseif checkLicense('bqpd')
      SolverFP='bqpd';
   elseif checkLicense('xa')
      SolverFP='xa';
   else
      SolverFP='milpSolve';
   end
   % SolverFP=GetSolver('fp',Prob.LargeScale);
   Prob.SolverFP = SolverFP;
else
   SolverFP=Prob.SolverFP;
end

% Setup Structure used in FP calls
ProbLP = CreateProbQP(Prob, 'fp', max(500,3*n), PriLev-3, optParam.wait);
ProbLP.QP.c = zeros(ProbLP.N,1);


%
%	Init
%

k = 0;	       % Iteration index. Initial value comp is step 0.
cL1Old=Inf;
Accept=0;

% sTrustr --- part 1


% Step 0: Initialization
gamma_1=0.1; gamma_2=0.5; gamma_3=2;
eta_1=0.01;  eta_2=0.25;  eta_3=0.75;
my_1=0.05;   my_2=0.1;

p_km1=zeros(n,1);
p=realmax*ones(n,1);

global n_f
if IterPrint
   fprintf('\n');
   fprintf('Iteration Function         f(x)           |step|     ');
   fprintf('trust region   |gradient|   Constraint\n');
   fprintf('             Count                                     ');
   fprintf('   radius       (projected)   Violations\n');
   gProj=realmax;
end

% sTrustr --- part 1 - end

set = zeros(n,1);

% Nonlinear constraints

if m > 0
   cL1 = norm(max(0,c_L-c_k),1) + norm(max(0,c_k-c_U),1);
else
   cL1=0;
   fMNew=Inf;
end

bl=cat(1,b_L,c_L);
bu=cat(1,b_U,c_U);
z =cat(1,Ax,c_k);

%iE=cat(1,xEqual,bEqual,cEqual);
%M=length(bl);
%I=find(~iE);
%E=find( iE);

u_k = -inf*ones(mTot,1);

% Check if x_0 is feasible

%h_k = norm(max(0,bl-z),1) + norm(max(0,z-bu),1);

zErr = max(0,max(bl-z,z-bu));
zSum = sum(zErr);

zErrOld = zErr;

if any(zErr > zTol .* max(1,abs(z)))
   if PriLev > 1
      xprinte(x_k,'x_k:');
      fprintf('Sum of infeasibilities             %30.20f\n',zSum);
   end

   ProbLP.A          = [A;dc_k];
   ProbLP.b_U        = bu-z;
   ProbLP.b_L        = bl-z;
   ProbLP.x_L        = x_L-x_k;
   ProbLP.x_U        = x_U-x_k;
   ProbLP.mLin       = size(ProbLP.A,1);

   ResultLP = tomRunMini(SolverFP,ProbLP);

   ExitFlag = ResultLP.ExitFlag;
   p        = ResultLP.x_k;
   v_k      = ResultLP.v_k;

   B        = ResultLP.QP.B;
   f_k      = ResultLP.f_k;


   if ExitFlag==0  % OK, feasible solution found.
      if PriLev > 2
         xprinte(p  ,'p:  ');
      end
      x_k = x_k + p;
      for i=1:M
          if M==1
             Prob.PartSep.index=0;
          else
             Prob.PartSep.index=i;
          end
          f_ik(i) = nlp_f(x_k, Prob, varargin{:});
      end
      f_k = sum(f_ik);
      if isinf(f_k)   % Too big step
         ExitFlag=1;
         x_k = x_k - p;
         for i=1:M
             if M==1
                Prob.PartSep.index=0;
             else
                Prob.PartSep.index=i;
             end
             f_ik(i) = nlp_f( x_k, Prob, varargin{:});
         end
         f_k = sum(f_ik);
      else
         if PriLev > 1
            fprintf('Phase I simplex found feasible point:\n');
            xprinte(x_k,'x_k:  ');
         end
         if mA > 0, Ax=A*x_k; end
         if ~isempty(c_k),
            c_k=nlp_c(x_k,Prob, varargin{:});  % Constraints
            dc_k = nlp_dc(x_k, Prob, varargin{:});
            cL1 = norm(max(0,c_L-c_k),1) + norm(max(0,c_k-c_U),1);
         end
         z=[Ax(:);c_k(:)];
         %h_k = norm(max(0,bl-z),1) + norm(max(0,z-bu),1);
         zErr = max(0,max(bl-z,z-bu));
         zSum = sum(zErr);
         if PriLev > 2
            fprintf('Sum of infeasibilities             %30.20f\n',zSum);
         end
      end
   end
end



cEqual=c_L==c_U;

stop = 0; f_km1=Inf;
alpha0cnt=0; alpha=1; fred0cnt=0; eig_step=0;

% sTrustr --- part 2

% Determine initial trust region radius

jMax=1;
iMax=5;
D=zeros(M,1);

x_0 = x_k;
for i = 1:M
    if M==1
       Prob.PartSep.index=0;
    else
       Prob.PartSep.index=i;
    end
    if alg==0
       [D(i),f_0i(i)] = itrr(x_k, alg, [],  jMax, iMax, Prob, varargin{:});
    else
       [D(i),f_0i(i)] = itrr(x_k, alg, B_k, jMax, iMax, Prob, varargin{:});
    end
end

f_0=sum(f_0i);
fM=f_0+zSum;
if m+mA > 0 & PriLev > 2
   fprintf('f_0 %40.25f\n',f_0);
   fprintf('fM    -     merit function value     %30.20f\n',fM);
end

Result.f_0 = f_0;

if PriLev >= 1 & any(x_k ~= x_0)
   disp('Start value updated by the initial trust region algorithm:')
   xprint(x_k,'x_k:');
   xprinte(D,'D:  ')
end

% Save the radius in D_vec
D_vec=D';

% Trust region constraint defined by inf-norm

f_ik=zeros(M,1);

% sTrustr --- part 2 END

while 1

   if k==0
      %for i=1:M
      %    Prob.PartSep.index=i;
      %    f_ik(i) = nlp_f( x_k, Prob, varargin{:});
      %end
      %f_k = sum(f_ik);

      f_ik=f_0i(:);
      f_k = f_0;

      Prob.PartSep.index=0;
      g_k = nlp_g(x_k, Prob, varargin{:});
      g_km1 = g_k;
   end
   if PriLev > 1
      fprintf('==================================================\n');
      fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);
      if mTot > 0
         fprintf('Sum of infeasibilities             %30.20f\n',zSum);
      end
   elseif IterPrint
      fprintf(' %5d  %7d   %20.17f %11.7e ',k, n_f, f_k, min(realmax,norm(p)));
      if max(abs(D)) < 1E-5
         fprintf('%13.5e',max(abs(D)));
      else
         fprintf('%13.5f',max(abs(D)));
      end
      fprintf('%16.7e %10.5e\n', norm(gProj), zSum);
   end

   if PriLev > 1
      xprint(x_k,'x_k:');
      xprint(set,'set:',' %2.0f',25);
      xprinte(g_k,'g_k:');
   end

%
% Compute the Hessian
%
   if alg==0 & k==0
      Prob.PartSep.index=0;
      H_k = nlp_H(x_k,Prob, varargin{:});
      if PriLev > 2 & alg==0 & m==0 & Prob.NumDiff==0
         fprintf('Eigenvalues at start\n');
         xprinte(eig(full(H_k)),'eig:');
      end
      if PriLev > 4 & alg==0
         if Prob.NumDiff==0
            fprintf('The Hessian H_k\n');
         else
            fprintf('The numerical Hessian approximation H_k\n');
         end
         PrintMatrix(H_k,'H_k:','%10.6e',6);
      end
   end

%
% Check convergence conditions
%

   Inform=0;

   % Check convergence criteria


   if max(abs(x_k-x_km1)./max(abs(x_k),size_x)) <= eps_x
      %& ...  m==0 | (m > 0 & max(D) < eps_x)
      if PriLev >= 1
         disp('*** Convergence 1, Iteration points are close ***');
      end
      Inform=Inform+1;
      flag=0;
      stop = 1;
   end
   Result.x_k=x_k;
   Result.g_k=g_k;
   Result.c_k=c_k;
   Result.cJac=dc_k;

   [v_k, Zv, P, Z_L, cErr, ceq, cineq, gProj] = LagMult(Prob,Result);

   Result.v_k=v_k;

   if m > 0 & PriLev > 2
      xprinte(cErr(1:min(length(c_k),MAX_c)),'cErr');
      if ceq > 0
         fprintf('%d equalities off more than cTol = %15.3e\n',ceq,cTol);
      end
      if cineq > 0
         fprintf('%d inequalities off more than cTol = %15.3e\n',...
                  cineq,cTol);
      end
   end

   m2=min(size(Z_L));
   if m2 < size(Z_L,2)
      fprintf('sTrustr found %d active constraints. ',size(Z_L,2))
      fprintf('Estimate Lagrange multipliers for %d first\n',m2);
   end
   if m2 > 0
      maxgProj = max(abs(gProj));
      % v = Z_L \ g_k;
      if PriLev > 2
         v = v_k(find(P~=0))
         fprintf('sTrustr Lagrange multiplier estimate for active constraints\n');
         xprinte(v,'v:');
      end
   else
      maxgProj = 0;
   end
   if any(maxgProj * max(abs(x_k),size_x) <= eps_g * max(abs(f_k),size_f)) ...
        & all(zErr <= zTol .* max(1,abs(z))) & ...
          all(max(D) < eps_x*max(1,abs(x_k)))

      if PriLev >= 1
         disp('*** Convergence 2, Projected gradient small ***');
      end
      Inform=Inform+2;
      flag=0;
      stop = 1;
   end

   if fred0cnt > LowIts      % Convergence criteria 3
      if PriLev >= 1
         fprintf('*** Convergence 3, Relative function value reduction')
         fprintf(' low for %d iterations ***\n',LowIts);
      end
      flag=0;
      Inform = Inform + 4;
      stop = 1;
   end
   if max(abs(D)) <= 0.001 * eps_x *  max([abs(x_k);size_x])
      if PriLev >= 1
         disp('*** Convergence 4, Too small trust region ***');
      end
      Inform=Inform+8;
      if any(zErr > zTol .* max(1,abs(z)))
         flag=1;
      else
         flag=0;
      end
      stop = 1;
   end
   % This test makes the algorithm stop too quick, before convergence
   %if max(abs(D)) > norm(p,inf) & all(zErr <= zTol .* max(1,abs(z))) & Accept
   %   if PriLev >= 1
   %      disp('*** Convergence 5, Can not proceed any longer ***');
   %   end
   %   Inform=Inform+16;
   %   flag=0;
   %   stop = 1;
   %end
   if k >= MaxIter		% Stop criteria 1
      if PriLev >= 1
         disp('*** STOP 1! Max no of iterations reached ***');
      end
      if Inform >0 & Inform < 100, flag=0;
      else
          flag=1; Inform=101;
      end
      stop = 1;
   end

   if f_k <= fLow       % Stop criteria 2
      if PriLev >= 1
         disp('*** STOP 2! Function value below given estimate   ***');
         disp('*** Restart with lower Prob.f_Low if min not reached ***');
      end
      if Inform >0 & Inform < 100, flag=0;
      else
          flag=1; Inform=102;
      end
      stop = 1;
   end

   if max(abs(D)) > 1E15
      if PriLev >= 1
         disp('*** STOP 3! Too large trust region, failure ***');
         disp('*** Problem non-convex, unbounded or ill-conditioned    ***');
      end
      if Inform >0 & Inform < 100, flag=0;
      else
          flag=1; Inform=103;
      end
      stop = 1;
   end


   if stop			% Show final results before return
      if PriLev >= 5  % REMARK! These results are displayed in the driver
         fprintf('==================================================\n');
         fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);

         xprint(x_k,'x_k:');
         xprinte(g_k,'g_k:');

         if ~isempty(H_k)
            if Prob.NumDiff ~= 0
               fprintf('The Numerical Hessian approximation H_k\n');
            else
               fprintf('The Hessian H_k\n');
            end
            mPrint(H_k,'H_k:',' %10.6e',6);
         end

      end
if DEBUG
   format long
   if mA > 0
      disp('[b_L A*x_k b_U]');
      [b_L A*x_k b_U]
   end
   if m > 0
      disp('[c_L c_k c_U]');
      [c_L c_k c_U]
   end
   format short
end
      Result.Iter=k;
      Result.ExitFlag=flag;
      Result.Inform=Inform;
      Result.ExitText=ExitText(Inform);
      Result.f_k=f_k;
      Result.B_k=B_k;
      Result.H_k=H_k;
      Result.x_0=x_0;
      Result.cJac=dc_k;
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=(x_k==x_L)+2*(x_k==x_U);
      Result=endSolve(Prob,Result);
      return;
   end
	
   %
   % Determine search direction p
   %
   % sTrustr --- Initial part in each iteration

   % Step 1: Model choice
   % (2.8) satisfied by definition

   % sTrustr --- Determining step by solving QPI

   % Step 2: Determination of step
   %% Solve for p_i, centered at x_k, initial trial p_km1

   if alg==0
      ProbQP.QP.F=H_k;
   else
      ProbQP.QP.F=B_k;
   end
   ProbQP.QP.c=g_k;
   ProbQP.x_L=max(-min(D)*ones(n,1),x_L-x_k);
   ProbQP.x_U=min(min(D)*ones(n,1),x_U-x_k);
%  ProbQP.x_0=max(ProbQP.x_L,min(ProbQP.x_U,p_km1));
   ProbQP.x_0=zeros(n,1);

   itMax = 5;
   it = 0;
   ExitFlag=2;
   while it < itMax+1 & ExitFlag~=0
   if mTot > 0
      ProbQP.A=[A;dc_k];
      ProbQP.b_L=bl-z;
      ProbQP.b_U=bu-z;
      if m > 0
         ix=find(c_k-c_U > cTol*max(1,abs(c_k)));
         % BAD FIX for nonlinear constraints
         ProbQP.b_U(mA+ix)=(1-itMax/100+it*0.01)*ProbQP.b_U(mA+ix);
         ix=find(c_k-c_L < cTol*max(1,abs(c_k)));
         % BAD FIX for nonlinear constraints
         ProbQP.b_L(mA+ix)=(1-itMax/100+it*0.01)*ProbQP.b_L(mA+ix);
      end
      ProbQP.mLin = size(ProbQP.A,1);
   else
      ProbQP.A=[];
      ProbQP.b_L=[];
      ProbQP.b_U=[];
   end

   ResultQP = tomRunMini(SolverQP,ProbQP);

   ExitFlag=ResultQP.ExitFlag;
   it = it + 1;
   end

   p=ResultQP.x_k;
   u_k=ResultQP.v_k;

if DEBUG
   disp('QP: x_L p x_U ; F  c ; A ; b_L A*p b_U')
   D
   [ProbQP.x_L p ProbQP.x_U]
   [ProbQP.QP.F ProbQP.QP.c]
   ProbQP.A
   [ProbQP.b_L ProbQP.A*p ProbQP.b_U]
   %pause
end
if DEBUG
   %x_k
   p
   u_k
   pause
end

   % sTrustr --- Determining step by solving QPI - END

   % sTrustr --- Calculate function changes

   m_k=f_k;
   if alg==0
      m_kp=f_k + g_k'*p + 0.5*(p'*H_k*p);
   else
      m_kp=f_k + g_k'*p + 0.5*(p'*B_k*p);
   end
   dm_k=m_k-m_kp;

   m_ik=f_ik;
   m_ikp=zeros(M,1);
   for i=1:M
       if M==1
          Prob.PartSep.index=0;
       else
          Prob.PartSep.index=i;
       end

       g_ik = nlp_g(x_k, Prob, varargin{:});

       if alg==0
          H_ik = nlp_H(x_k, Prob, varargin{:});
       else
          H_ik = B_k;
       end

       m_ikp(i) = f_ik(i) + g_ik'*p +  0.5*(p'*H_ik*p);
   end
   dm_ik=m_ik-m_ikp;

   f_ikm1=f_ik;

   % New function value in x_k + p
   x_kp=x_k + p;
   if PriLev > 1
      xprinte(p,'p_?:');
      xprinte(x_kp,'x_?:');
   end

   for i=1:M
       if M==1
          Prob.PartSep.index=0;
       else
          Prob.PartSep.index=i;
       end
       f_ikp(i) = nlp_f(x_kp, Prob, varargin{:});
   end
   f_ikp=f_ikp(:);
   f_kp = sum(f_ikp);

   df_k=f_k-f_kp;

   fp_0=g_k'*p;
   if PriLev > 2
      fprintf('Directed derivative: %15.7e\n',fp_0);
   end

   if m > 0
      c_k  = nlp_c( x_kp, Prob, varargin{:});
      cL1 = norm(max(0,c_L-c_k),1) + norm(max(0,c_k-c_U),1);
      fMNew=f_kp+cL1;
      if PriLev > 2
         fprintf('cL1 - nonlinear constraint violation %30.20f\n',cL1);
         fprintf('fMNew - new merit function value     %30.20f\n',fMNew);
      end
   end

   % Step 3: Overall model fit. Should the step p be accepted?
   if (cL1 <=cTol & df_k >= eta_1*dm_k ...
      | (fred0cnt > 3 & ~isinf(df_k) & df_k > 0))...
      | fMNew <= (1-1E-10)*fM  | ...
      ( m > 0 & M==1 & M(1) > norm(p,inf) & cL1 <= cL1Old )
      % SHIT| (m > 0 & (zSum < zSumOld) | df_k > 0)
      % (2.43) + check on merit function

      if PriLev >=1
         if f_kp > f_k & m == 0
            disp('sTrustr: Function value increases !!! No descent?)')
         end
      end
      Accept=1;
      f_km1=f_k;
      x_km1=x_k;
      g_km1=g_k;
      x_k=x_kp;                % (2.44)
      p_km1=p;
      f_k=f_kp;
      f_ik=f_ikp;
      Prob.PartSep.index=0;
      g_k  = nlp_g(x_k, Prob, varargin{:});
      if alg==0
         H_k  = nlp_H(x_k, Prob, varargin{:});
      end
      if m > 0
         dc_k = nlp_dc(x_k, Prob, varargin{:});
      end
      if mA > 0, Ax=A*x_k; end
      if mTot > 0
         z=[Ax(:);c_k(:)];
         %h_k = norm(max(0,bl-z),1) + norm(max(0,z-bu),1);
         zErr = max(0,max(bl-z,z-bu));
         zSum = sum(zErr);
      end
      % J=zeros(mTot,1);
      % Check on constraint violation
      % J(bu + cTol < z)= 1;              % Upper bounds does not hold
      % J(bl - cTol > z)=-1;              % Lower bounds does not hold
   else
      Accept=0;
      if PriLev > 1
         disp('DO NOT Accept step');
      end
   end

   df_ik=f_ikm1-f_ikp;

   % Step 4: Update element trust region radii
   for i=1:M
      if abs(dm_ik(i)) <= my_1*dm_k/M % Case 2: negligible element
         if PriLev > 2
            fprintf('Negligible. Element %d\n',i)
            element=i;
         end
         if (abs(df_ik(i)) <= my_2*dm_k/M) & (df_k >= eta_1*dm_k)
%disp('C1')
            D(i)=gamma_3*D(i);      % (2.57)
         elseif abs(df_ik(i)) <= my_2*dm_k/M
%disp('C2')
            % Do nothing                      % (2.58)
         else
%disp('C3')
            %choice in [gamma_1*D(i),gamma_2*D(i)];
            D(i)=gamma_1*D(i);
         end
      else                       % Case 1: meaningful element
         if PriLev > 2
            fprintf('Meaningful. Element %d\n',i)
            element=i;
         end
         if (df_ik(i) >= dm_ik(i)-(1-eta_3)*dm_k/M) & (df_k >= eta_1*dm_k)
%disp('D1')
            D(i)=gamma_3*D(i);      % (2.51)
         elseif (df_ik(i) >= dm_ik(i)-(1-eta_3)*dm_k/M) | it > 1 ...
            % Do nothing                      % (2.52)
            if ~Accept & m > 0 & M==1
               % Must do something!
               D(i)=gamma_1*D(i);      % (2.51)
            end
         elseif df_ik(i) >= dm_ik(i)-(1-eta_2)*dm_k/M
            D(i)=gamma_2*D(i);      % (2.54)
%disp('D3')
         else
%disp('D4')
            D(i)=gamma_1*D(i);
            %if all( zErr <= zTol .* max(1,abs(z))) | fMNew > (1-1E-1)*fM
%format long
%disp('fm-fmnew')
%fM-fMNew
%pause
            %   if fM-fMNew < -1E-7 & all(zErr <= zTol .* max(1,abs(z)))
            %      % Here we have trouble
            %      D(i)=gamma_1*D(i);
            %   else
            %      D(i)=gamma_1*D(i);
            %   end
            %elseif m>0 & fMNew < 0.1*fM
            %   %disp('EXPAND');
            %   D(i)=gamma_3*D(i);      % (2.57)
            %end
         end
      end
   end

if PriLev > 1
   xprinte(D,'D:  ')
end

   % sTrustr --- END

   if m == 0 & df_k <= eps_f * max(f_k,size_f);
      fred0cnt=fred0cnt+1;
   else
      fred0cnt=0;
   end

   if optParam.wait, pause; end;
%
% Update
%

   if Accept
      fM=fMNew;
      zErrOld=zErr;
      if alg > 0  % Only update Quasi-Newton if step is accepted
         B_k=SafeUpdate(alg,alpha,p,xTol,idx,g_k,g_km1,B_k,PriLev,k);
      end
   end

   cL1Old=cL1;
   k = k+1;

end	% (while)



% =========================================
function [alg, SolvAlg, B_k, LowIts]=NameAlg(Prob,M,userH,PriLev)
% =========================================

% alg = 0 = Structured Trust Region using Hessian, numerical or analytical
% alg = 1 = Structured Trust Region, BFGS update
% alg = 2 = Structured Trust Region, DFP update

alg=max(0,Prob.Solver.Alg);
if isempty(alg), alg=0; end
if alg > 2, alg=0; end

if alg==0 & Prob.NumDiff==0 & (~isempty(userH) | any(Prob.ADObj == -1))
   SolvAlg='Structured Trust Region Algorithm';
elseif alg==0
   SolvAlg='Structured Trust Region with Numerical Hessian';
elseif alg==1
   SolvAlg='Structured Trust Region with BFGS update';
elseif alg==2
   SolvAlg='Structured Trust Region with DFP  update';
end

if M > 1
   SolvAlg = [SolvAlg '; ' num2str(M) ' psf'];
end

% No of iterations with low reduction before convergence
LowIts = 10;

% Find initial Quasi-Newton matrix
n = length(Prob.x_0);
if alg~=0
   if ~isempty(Prob.optParam.QN_InitMatrix)
      B_k = Prob.optParam.QN_InitMatrix;   % B = User given matrix
      if size(B_k,1)~=n | size(B_k,2)~=n
         if PriLev >=0
            disp('Illegal dimensions on initial Quasi Newton matrix');
            disp('sTrustr is using identity matrix instead');
         end

         B_k = eye(n);                % B = Use identity matrix
      end
   else
      B_k = eye(n); % B = Start Quasi-Newton with identity matrix
   end
   LowIts = 30; % Algorithm is slow for QN. Do not stop on low reduction now
else
   B_k=[];
end

% =========================================
function [D_0, f_0, x_0] = itrr(x_0, alg, B_k, jMax, iMax, Prob, varargin)
% =========================================

%
% itrr.m = ITRR - Initial Trust Region Radius
% Utility used by strustr.m for Structured Trust Region Algorithm
% Ref. A. Sartenaer, Report 95/4, June 16, 1995
%
%
% function [D_0, f_0, x_0] = itrr(x_0, alg, B_k, jMax, iMax, Prob, varargin)
%
% INPUT PARAMETERS
%
%    x_0    Starting point
%    alg    Algorithm
%    B_k    Hessian approximation if alg > 0
%    jMax   Number of outer iterations, normally 1
%    iMax   Number of inner iterations, normally 5
%    Prob   TOMLAB problem structure
%
%    Prob.PartSep.index Index for the partial function to be analyzed
%
%    varargin - Extra user parameters, passed to f,g,and H; e.g. user_P
%
% OUTPUT PARAMETER
%    D_0    = Initial trust region radius
%    f_0    = Function value at the input point x_0
%    x_0    = Updated starting point, if jMax > 1
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written Sept 12, 1996. Last modified June 26, 1999.
%


if isempty(jMax), jMax=1; end
if isempty(iMax), iMax=5; end

% Step 0: Initialization of starting point loop

gamma_1=0.0625; gamma_2=5; gamma_3=0.5; gamma_4=2;

my_0=0.5; my_1=0.5;my_2=0.35;

teta=0.25;

% Step 1: Starting point loop

for j=0:jMax
    x_k = x_0;
    f=nlp_f(x_k, Prob, varargin{:});

    if j==0, f_0=f; end    % Save initial f value as f_0

    g=nlp_g(x_k, Prob, varargin{:});

    if norm(g,2) > eps
       g_n=g/norm(g,2);
    else
       g_n=g;
    end
    if alg==0
       B=nlp_H(x_k, Prob, varargin{:});
    else
       B=B_k;
    end
    sigma=0;
    delta=0;

    % Step 1.0: Initialization of radius loop
    % itr=norm(0.1*x_0,2);
    itr=1;
    itrmax=0;

    % Step 1.1: Radius update loop

    for i=0:iMax
        s=itr*g_n;
        m_i=f-g'*s-0.5*(s'*B*s);

        x_k=x_0-s;
        f_i=nlp_f(x_k, Prob, varargin{:});

        if f-m_i == 0  % safeguard against division by 0
           rho=Inf;
        else
           rho=(f-f_i)/(f-m_i);
        end
        rho1=abs(rho-1);

        % Step 1.2: Maximum radius estimate update
        if rho1 <= my_0 itrmax=max(itrmax,itr); end;

        % Step 1.3: Best decrease update
        if j < jMax &  f-f_i > delta
           delta=f-f_i;
           sigma=itr;
        end
        % Step 1.4 Radius estimate update
        if i < iMax
           b1=(teta*(f-g'*s)+(1-teta)*m_i-f_i);
           if b1==0
              beta(1)=Inf;
           else
              beta(1)=-teta*g'*s/b1;    % (2.11)
           end
           b2=(-teta*(f-g'*s)+(1+teta)*m_i-f_i);
           if b2==0
              beta(2)=Inf;
           else
              beta(2)=teta*g'*s/b2;    % (2.12)
           end
           if rho1 > my_1
              if min(beta) > 1                                     % (2.13)
                 beta_i=gamma_3;
              elseif max(beta) < gamma_1 | (min(beta) < gamma_1 & max(beta) > 1)
                 beta_i=gamma_1;
              elseif (beta(1) >= gamma_1 & beta(1) <= 1) & ...
                     (beta(2) <  gamma_1 | beta(2) >  1)
                 beta_i=beta(1);
              elseif (beta(2) >= gamma_1 & beta(2) <= 1) & ...
                     (beta(1) <  gamma_1 | beta(1) >  1)
                 beta_i=beta(2);
              else
                 beta_i=max(beta);
              end
           elseif rho1 <= my_2
              if max(beta) < 1                                     % (2.14)
                 beta_i=gamma_4;
              elseif max(beta) > gamma_2
                 beta_i=gamma_2;
              elseif (beta(1) >= 1 & beta(1) <= gamma_2) & (beta(2) < 1)
                 beta_i=beta(1);
              elseif (beta(2) >= 1 & beta(2) <= gamma_2) & (beta(1) < 1)
                 beta_i=beta(2);
              else
                 beta_i=max(beta);
              end
           else
              if max(beta) < gamma_3                                % (2.15)
                 beta_i=gamma_3;
              elseif max(beta) > gamma_4
                 beta_i=gamma_4;
              else
                 beta_i=max(beta);
              end
           end
        end
    end

    % Step 2: Starting point update
    if j < jMax & delta > 0
       x_0=x_0-sigma*g_n;
    end
end

% Step 3: Final radius

if itrmax > 0
   D_0=itrmax;
else
   D_0=itr;
end

% MODIFICATION LOG:
%
% 981020  hkh   Changed to use Prob structure for sending info to functions
% 981021  hkh   Safeguard rho computation against division by 0, f-m_i==0
% 981028  hkh   Compute f_0 = f(x_0) and return as 2nd argument
% 990626  hkh   Avoid feval. Accept numerical approximation of Hessian

% =========================================
function B_k=SafeUpdate(alg,alpha,p,xTol,idx,g_k,g_km1,B_k,PriLev,k)
% =========================================
   if alg == 1 % Safeguarded BFGS update of approximate Hessian
      z=p;
      if norm(z) > xTol
         y=g_k(idx)-g_km1(idx);
         Bz=B_k(idx,idx)*z;
         zBz=z'*Bz;
         zy=z'*y;
         if zy < 0.2 * zBz
            theta=0.8*zBz / (zBz - zy);
            w=theta * y + (1-theta)*Bz;
            zw=z'*w;
         else
            w=y;
            zw=zy;
         end
         if PriLev > 2
            fprintf('BFGS Hessian update iter %3.0f\n',k);
            fprintf('zw %10.6e zBz %10.6e \n',zw,zBz)
         end
         if zw < 1E-13 | zBz < 1E-13
            if PriLev > 1
               fprintf('BFGS Hessian update step dangerous!\n');
               fprintf('max(w)/zw %20.10e\n',max(w)/zw);
               fprintf('min(w)/zw %20.10e\n',min(w)/zw);
               fprintf('max(Bz)/zBz %20.10e\n',max(Bz)/zBz);
               fprintf('min(Bz)/zBz %20.10e\n',min(Bz)/zBz);
            end
         end
         if zw == 0
            if zBz ~=0
               B_k(idx,idx)=B_k(idx,idx)-(Bz/zBz)*Bz';
            end
         elseif  zBz == 0
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w';
         else
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w'-(Bz/zBz)*Bz';
         end
         if PriLev > 3
            fprintf('B_k matrix after BFGS update\n');
            PrintMatrix(B_k,'B_k:')
            fprintf('B_k matrix eigenvalues:\n');
            xprint(eig(B_k),'eig:');
         end
      end
   elseif alg == 2 % Safeguarded DFP update of the approximate Hessian
      z=p;
      if norm(z) > xTol
         y=g_k(idx)-g_km1(idx);
         Bz=B_k(idx,idx)*z;
         zy=z'*y;
         if PriLev > 2
            fprintf('DFP Hessian update iter %3.0f\n',k);
            fprintf('zy %10.6e\n',zy)
         end
         if zy < 1E-13
            if PriLev > 1
               fprintf('DFP Hessian update step dangerous!\n');
               fprintf('max(y)/zy %20.10e\n',max(y)/zy);
               fprintf('min(y)/zy %20.10e\n',min(y)/zy);
            end
         end
         if zy ~= 0
            B_k(idx,idx)=B_k(idx,idx)+((1+(z/zy)'*Bz)/zy)*y*y'-(y*Bz'+Bz*y')/zy;
         end
         if PriLev > 3
            fprintf('B_k matrix after BFGS update\n');
            PrintMatrix(B_k,'B_k:')
            fprintf('B_k matrix eigenvalues:\n');
            xprint(eig(B_k),'eig:');
         end
      end
   end

% ------------------------------
function Text = ExitText(Inform)
% ------------------------------

switch  Inform
   case 1
     Text = 'Iteration points are close';
   case 2
     Text = 'Projected gradient small';
   case 3
     Text = 'Iterations close, projected gradient small';
   case 4
     Text = 'Relative f(x) reduction low for LowIts iterations';
   case 5
     Text = 'Iterations close, relative f(x) reduction low';
   case 6
     Text = 'Proj grad small, f(x) reduction low for LowIts iterations';
   case 7
     Text = str2mat('Iterations close, projected gradient small' ...
                   ,'f(x) reduction low for LowIts iterations');
   case 8
     Text = 'Trust region small';
   case 9
     Text = 'Trust region small. Iteration points close';
   case 10
     Text = 'Trust region and projected gradient small';
   case 11
     Text = 'Trust region and projected gradient small, iterations close';
   case 12
     Text = 'Trust region small, Relative f(x) reduction low';
   case 13
     Text = str2mat('Trust region small, Relative f(x) reduction low' ...
                   ,'Iteration points are close');
   case 14
     Text = str2mat('Trust region small, Relative f(x) reduction low' ...
                   ,'Projected gradient small');
   case 15
     Text = str2mat('Trust region small, Relative f(x) reduction low' ...
                   ,'Iteration points close, Projected gradient small');
   case 101
     Text = 'Maximal number of iterations reached';
   case 102
     Text = str2mat('Function value below given estimate' ...
                   ,'Restart with lower fLow if minimum not reached');
   case 103
     Text = str2mat('Too large trust region' ...
                   ,'Problem non-convex, unbounded or ill-conditioned');
end
if Inform < 100
   Text=str2mat('Optimal solution found',Text);
end

% MODIFICATION LOG:
%
% 981021  hkh   Revised code, using ucSolve, for new TOMLAB design
% 981023  hkh   Bug fixes. Deleting some ucSolve code, not relevant
%               Made it work for linear constraints.
% 981026  hkh   Check empty A, Add Ax(:), otherwise crash when empty matrix
% 981026  hkh   Solver is name of solver. SolverAlgorithm is description
%               Fix printing and algorithm choices
% 981027  hkh   Use Prob.Solver.Alg to define alg
% 981028  hkh   Put f_0 = f(x_0) as field in Result
% 981102  hkh   Add handling of nonlinear constraints.
% 981105  hkh   Unnecessary to set Result.Prob here, done in endSolve
% 981107  hkh   solvType=1 should be solvType=3, now when handling constraints
% 981110  hkh   Result.cJac instead of Result.dc_k. Use MAX_c and MAX_x
% 981120  hkh   Add FP to find feasible x_0.
% 981120        Correct the setting of x_0 before QP.
% 981124  hkh   Avoid bad starting value with infinite value
% 981129  hkh   Change printing levels for many printings
% 990820  hkh   Revision for 2.0.
% 990830  hkh   Use general relative tolerances. Check each constraint.
% 000916  hkh   Add ExitText information
% 000927  hkh   Must check on large trust region to avoid loop to MaxIter
% 001021  hkh   Add #psf to algorithm text string
% 001101  hkh   Use qld as default QP solver for v2.1 and v3.0 /MINI
% 031201  hkh   Revise AD, new field ADObj, Revise algorithm selection
% 040111  hkh   Change call to inisolve
% 040125  hkh   Define field mLin in ProbLP and ProbQP
% 050117  med   mlint revision
% 051216  med   Help updated with Inform values and text
% 060814  med   FUNCS used for callbacks instead
% 060818  hkh   Use Prob.f_Low instead of optParam.eps_absf for convergence test
% 061212  med   ADMAT removed
% 070907  hkh   SolverQP/FP picked from list, 1st with license, avoid GetSolver
% 080607  hkh   Use tomRun, not tomSolve
% 080607  med   Switched to tomRunMini
% 090717  med   QP calculations updated
