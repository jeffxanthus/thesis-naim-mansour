% clsSolve implements algorithms for weighted nonlinear least
% squares problems with linear inequality and equality constraints and
% bound constraints.
%
% clsSolve solves the (weighted) nonlinear least squares problem of the form:
%
%	min       0.5 * sum {r(x).^2} <==>      min       0.5 *  r(x)^T * r(x)
%        x                                       x
%                 x_L <=   x  <= x_U                      x_L <=   x  <= x_U
%                 b_L <=  Ax  <= b_U                      b_L <=  Ax  <= b_U
%
% where r is the m-dimensional residual vector and
% x is a n-dimensional unknown parameter vector
%
% Residual weighting may be applied in a general way, see help clsAssign
% and the parameters weightType and weightY.
%
% function Result = clsSolve(Prob, varargin)
%
% Algorithms:
% Prob.Solver.Alg = 0, empty or illegal number, gives the Fletcher-Xu algorithm
%
% Prob.Solver.Alg = 1  The Fletcher-Xu hybrid method; Gauss-Newton / BFGS
% Prob.Solver.Alg = 2  The Al-Baali - Fletcher hybrid method; Gauss-Newton/BFGS
% Prob.Solver.Alg = 3  Huschens method. SIAM J. Optimization. Vol 4, No 1,
%                      pp 108-129 jan 1994.
% Prob.Solver.Alg = 4  The Gauss-Newton method with subspace minimization
% Prob.Solver.Alg = 5  Wang, Li, Qi Structured MBFGS method
% Prob.Solver.Alg = 6  Li-Fukushima MBFGS method
% Prob.Solver.Alg = 7  Broydens method
%
% Recommendations: Alg=5 is theoretically best, and seems best in practice as
% well. Alg=1 and Alg=2 behave very similar, and are robust methods.
% Alg=4 may be good for ill-conditioned problems.
% Alg=3 and Alg=6 may sometimes fail
% Alg=7 tries to minimize Jacobian evaluations, but might need more
% residual evaluations. Also fails more often that other algorithms
% Suitable when analytic Jacobian is missing and evaluations of the Jacobian
% is costly. The problem should not be too ill-conditioned
%
% Search method technique (if Prob.LargeScale = 0):
% Prob.Solver.Method = 0 QR with pivoting (both sparse and dense)
% Prob.Solver.Method = 1 SVD (dense)
% Prob.Solver.Method = 2 Matlabs inversion routine (Uses QR) (unreliable)
% Prob.Solver.Method = 3 Explicit computation of pseudoinverse, pinv(J_k) (unreliable)
%
% Search method technique (if Prob.LargeScale = 1, then Method = 0 always):
% Prob.Solver.Method = 0 Sparse iterative QR using Tlsqr
%
% Bound constraints treated as described in
% Gill, Murray, Wright: Practical Optimization, Academic Press, 1981.
% Null space method for the equality constraints.
%
% clsAssign.m may be used to define the Prob structure
%
% INPUT PARAMETERS
%
% Prob    Structure, where the following variables are used:
%
%   Solver.Alg     Described above
%   Solver.Method  Described above
%
%   LargeScale If = 1, then sparse iterative QR using Tlsqr is used to find
%              search directions
%
%   x_0       Starting point
%   x_L       Lower bounds for x
%   x_U       Upper bounds for x
%   b_L       Lower bounds for linear constraints
%   b_U       Upper bounds for linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for nonlinear constraints
%   c_U       Upper bounds for nonlinear constraints
%   f_Low     A lower bound on the optimal function value,
%             see LineParam.fLowBnd below
%
%   SolverQP  Name of the solver used for QP subproblems. If empty,
%             picked from a list, best available with a license
%
%   PriLevOpt Print Level.
%   optParam Optimization parameters. Fields used in optParam are:
%      bTol, eps_absf, eps_g, eps_Rank, eps_x, IterPrint, MaxIter,
%      PreSolve size_f, size_x, xTol, wait
%
%      optParam.QN_InitMatrix:  Initial Quasi-Newton matrix, if not empty,
%                                 otherwise use identity matrix
%   LineParam Line search parameters, special fields used:
%
%   if alg == 7:
%      LineAlg  = 0 Fletcher quadratic interpolation line search
%      LineAlg  = 3 Fletcher cubic interpolation line search
%      otherwise Armijo-Goldstein line search (LineAlg == 2)
%   if alg ~= 7:
%      LineAlg  = 0 Fletcher quadratic interpolation line search
%      LineAlg  = 1 Fletcher cubic interpolation line search
%      LineAlg  = 2 Armijo-Goldstein line search
%      otherwise  Fletcher quadratic interpolation line search (LineAlg == 0)
%
%   If Fletcher, see help LineSearch for the LineParam parameters used
%   Most important is the accuracy in the line search
%   sigma        Line search accuracy tolerance, default 0.9
%
%   If Armijo-Goldstein, then 2 parameter fields in LineParam is used:
%   agFac        Armijo Goldsten reduction factor, default 0.1
%   sigma        Line search accuracy tolerance, default 0.1
%
%   fLowBnd      A lower bound on the global optimum of f(x).
%                NLLS problems always have f(x) values >= 0
%                The user might also give lower bound estimate in Prob.f_Low
%                clsSolve computes LineParam.fLowBnd as:
%                LineParam.fLowBnd = max(0,Prob.f_Low,Prob.LineParam.fLowBnd)
%                fLow = LineParam.fLowBnd is used in convergence tests
%
% Extra parameters:
% VARARGIN: User defined parameters passed to user defined functions
%
% OUTPUT PARAMETERS
%
% Result      Structure with results from optimization
%    x_k      Optimal point
%    v_k      Lagrange multipliers NOT USED
%    f_k      Function value at optimum
%    g_k      Gradient vector at optimum
%    x_0      Starting value vector
%    f_0      Function value at start
%    r_k      Residual at optimum
%    J_k      Jacobian matrix at optimum
%    xState   Variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%    bState   Linear constraint: Inactive==0; On lower bound == 1;
%             On upper bound == 2; Equality == 3;
%    Iter     Number of iterations
%    ExitFlag Flag giving exit status
%       0     Convergence
%       1     Error, see Inform
%    ExitTest Text string giving ExitFlag and Inform information
%    Inform   Code telling type of convergence
%             1   Iteration points are close.
%             2   Projected gradient small.
%             3   Iteration points are close and projected gradient small.
%             4   Function value close to 0.
%             5   Iteration points are close and function value close to 0.
%             6   Projected gradient small and function value close to 0.
%             7   Iteration points are close, projected gradient small and
%                 function value close to 0.
%             8   Relative function value reduction low for 10 iterations.
%            11   Relative f(x) reduction low for LowIts iter. Close Iters.
%            16   Small Relative f(x) reduction.
%            17   Close iteration points, Small relative f(x) reduction.
%            18   Small gradient, Small relative f(x) reduction.
%
%            32   Local min with all variables on bounds.
%            99   The residual is independent of x. The Jacobian is 0.
%           101   Max no of iterations reached.
%           102   Function value below given estimate.
%           104   x_k not feasible, constraint violated.
%           105   The residual is empty. No NLLS problem
%
%    Solver   Solver used.
%    SolverAlgorithm    Solver algorithm used.
%    Prob     Problem structure used.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Jan 15, 1998.   Last modified Jul 22, 2011.

function Result = clsSolve(Prob, varargin)

%#function ls_f ls_g ls_H lls_r lls_J lls_H

if nargin < 1
   error('clsSolve needs input structure Prob');
end

% NOT USED YET:
%   FUNCS.d2r:The routine to compute the 2nd derivative wrt r, as a string
%   FUNCS.c:  The routine to evaluate the constraints as a string
%   FUNCS.dc: The routine to compute the gradient of the constraints

solvType=checkType('cls');

Prob=ProbCheck(Prob,'clsSolve',solvType);

Prob = iniSolve(Prob,solvType,1,1);

% If robust==1, then use more safe convergence tests
Robust=0;

DEBUG=0;     % Debugging flag

% If true (==1) plot line search problem at each iteration
plotLine=Prob.plotLine;

optParam  = Prob.optParam;
LineParam = Prob.LineParam;


f_red  = Inf;

[alg, SolvAlg, B_k, A_k]=NameAlg(Prob);

% HKH Add LineAlg parameter
if alg==7
   % Default Armijo for Broydens to avoid Jacobian evaluations
   LineAlg  = DefPar(Prob.LineParam,'LineAlg',1);
   if LineAlg == 3
      SolvAlg=[SolvAlg '. Fletcher cubic line search'];
      LineAlg = 1;
   elseif LineAlg == 0
      SolvAlg=[SolvAlg '. Fletcher quadratic line search'];
      LineAlg = 0;
   else
      SolvAlg=[SolvAlg '. Armijo non-derivative line search'];
      LineAlg = 2;
   end
   Prob.LineParam.LineAlg = LineAlg;
   LineParam.LineAlg = LineAlg;
   NewJ = 0;
else
   LineAlg  = DefPar(Prob.LineParam,'LineAlg',1);
   if LineAlg > 2, LineAlg = 1; end
end
% HKH agFac and sigma are used for Armijo-Goldstein line search
if LineAlg > 1
   agFac  = DefPar(LineParam,'agFac',0.1);
   sigma  = DefPar(LineParam,'sigma',0.1);
else
   sigma  = DefPar(LineParam,'sigma',0.9);
end

Result=ResultDef(Prob);
Result.Solver='clsSolve';
Result.SolverAlgorithm=SolvAlg;


[n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

p_rank =n;
pRank  =n;
pCond  =NaN;

if optParam.PreSolve & ~isempty(Prob.A)
   Prob     = preSolve(Prob);  % Call presolve analysis routine preSolve.m
   NaNidx   = find(~isnan(Prob.b_L)); % Constraints that can not be deleted.
   Prob.A   = Prob.A(NaNidx,:);
   Prob.b_L = Prob.b_L(NaNidx);
   Prob.b_U = Prob.b_U(NaNidx);
   Prob.mLin = size(Prob.A,1);
end

% Linear constraints

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

% Nonlinear constraints

[m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin{:});

if m > 0
   if ~all(c_L==c_U)
      disp('clsSolve do not handle nonlinear inequality constraints');
      disp('Run NLSSOL in TOMLAB /NPSOL (or TOMLAB /SOL)');
      disp('Or run slsSolve (using subsolver SNOPT in TOMLAB /SNOPT)');
      Result.ExitFlag = 1;
      Result.ExitText = 'Cannot solve problems with nonlinear inequality constraints';
      Result=endSolve(Prob,Result);
      return;
   else
      disp('clsSolve do not handle nonlinear equality constraints');
      disp('Try adding them as additional residuals, with a weight factor');
      disp('Or run NLSSOL in TOMLAB /NPSOL (or TOMLAB /SOL)');
      disp('Or run slsSolve (using subsolver SNOPT in TOMLAB /SNOPT)');
      Result.ExitFlag = 1;
      Result.ExitText = 'Cannot solve problems with nonlinear equality constraints';
      Result=endSolve(Prob,Result);
      return;
   end
end

mTot=mA+m;

lam_L = zeros(n,1);  % Initial zeros for Lagrange multipliers
%u_k = -inf*ones(mTot,1);

bl=cat(1,b_L,c_L);
bu=cat(1,b_U,c_U);
z =cat(1,Ax,c_k(:));

zEqual=cat(1,bEqual,cEqual);

zErr = max(0,max(bl-z,z-bu));
zSum = sum(zErr);

% Pick up input variables from Prob structure
Ud2r = Prob.FUNCS.d2r;

% Set to 1 if just releasing one variable at a time.
NOT_release_all = 0;

bTol    = optParam.bTol;     % Linear constraint feasibility tolerance
epsRank = optParam.eps_Rank; % Rank test tolerance
eps_x   = optParam.eps_x;    % Convergence tolerance in x
eps_absf = optParam.eps_absf;% Convergence tolerance in f close to 0
MaxIter = optParam.MaxIter;  % Maximal number of iterations
size_x  = optParam.size_x;   % Approximate size of optimal variable values
size_f  = optParam.size_f;   % Approximate size of optimal function value

% Lower bound on function value, not critical that it is tight
fLow              = max([0,Prob.LineParam.fLowBnd,Prob.f_Low]);
LineParam.fLowBnd = fLow;

% If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bounds
xTol =  optParam.xTol;
% If (f_km1-f_k) < eps_f * max(f_k,size_f) for LowIts == 10 iter, stop
%NOT USED NOW eps_f = optParam.eps_f;
% No of iterations with low reduction before convergence
LowIts = 10;

%ONLY USED if active set handling of NONLINEAR CONSTRAINTS: cTol=optParam.cTol;
%zTol=[bTol*ones(mA,1);cTol*ones(m,1)];

zTol=bTol*ones(mA,1);

%LineParam.eps1=min(1E-7,0.1*eps_x);
% WHY IS THIS BETTER ??? (Handles ill-scaled problems ???)
LineParam.eps1=1E-11;

alphaMax = 1E20;       % Set here for conv. check in first iteration.

Method    = Prob.Solver.Method;   % Search direction solution technique
PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration

if MaxIter <= 0, MaxIter=100*n; end

[Prob]=CheckSep(Prob);

%if 0 & LineParam.LineAlg > 1
%   CURV=1;
%   LineParam.LineAlg=LineParam.LineAlg-2;
%   global xLast
%   xLast=x_km1;
%else
   CURV=0;   % No curvi linear search
%end

LargeScale = Prob.LargeScale;

if isempty(Prob.SolverQP)
   % Convex = 1;
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
   % SolverQP=GetSolver('qp',LargeScale,Convex);
   Prob.SolverQP = SolverQP;
else
   SolverQP=Prob.SolverQP;
end

if n < 10
   LargeScale = 0;
end
% Always use Method 0 if LargeScale true
if LargeScale > 0, Method = 0; end

TESTONX=1;
KT=Inf;
KTT=0;

k = 0; % Iteration index
set = zeros(n,1);

xR=[]; xRp=[];
nract = 0; % Number of active variables
nractcon = 0; % Number of active inequality constraints


[x_k, me, mc, mi, LCflag, A, b, W_kset]=leqCheck(Prob, ...
      SolverQP, x_k, b_L, b_U, x_L, x_U, bEqual, Ax, bTol, PriLev);

if ~LCflag, Z=[]; end
Arank = me; % Pseudo rank for the active constraints. Assume full rank at x_0

if mc > 0 % If any linear constraints

   W_km1set = W_kset;
   W_k = find(W_kset); % Set of active constraints including all equalities
   if mi > 0
      actIcset = zeros(mc,1);
      actIc = find(actIcset); % Set of active inequality constraints
   end
   %Aw_k = A(W_k,:);
else
   W_kset = [];
end

lambda = 0; % Initial zeros for Lagrange multipliers
pNorm  = 0;
relcon = -1;
actcon = -1;
actvar = -1;
releaseall = 1;
xRR = []; % Variables to be released if releaseall
xRRlogic = [];

%r_k = nlp_r(x_k,Prob, varargin{:});   % Function value
%J_k = nlp_J(x_k,Prob, varargin{:});   % Jacobian
%f_k = 0.5 * (r_k' * r_k);
%g_k = J_k' * r_k;
%JJ=  J_k'*J_k;
%m = length(r_k);

H_k= [];
H2 = [];


IterGN = 0;
GN=1; % In Fletcher-Xu hybrid (alg==1) and Al-Baali-Fletcher hybrid (alg==2)
      % GN=1 ==> Take GN-step. GN=0 ==> Take BFGS-step.

alpha0cnt=0; stop = 0; alpha=1; f_km1=Inf; fRed0cnt=0;

% Initial dummy value
fp_0=0;
if IterPrint
   fprintf('Iteration Function         f(x)           |step|     ');
   fprintf('line search   |gradient|\n');
   fprintf('           Count                                       ');
   fprintf('   step        (Projected)\n');
end
p_full=0;
gProj=0;
alpha=0;

global n_r
while 1
% Check if variables near bound and set up active variable indicator set

   lambda = 1;

   %nract = 0;
   changed = 0;
   %x_z = x_k;
   % Adjust x_k inside bounds
   x_k = max(x_L,min(x_U,x_k));
   %changed = any(x_k~=x_z);

   if (mi <= 0) | k==0
      set_0=set;
      [set, nract]=x_kActive(nract, x_k, x_L, x_U, xTol, Arank);
   end

   if k==0,set_0=set;end % 1st iteration release from bounds directly
   b_idx = find(set);    % Index vector for active variables, the fixed vars
   idx = find(~set);     % Index vector for non active variables, the free vars
   if changed | k==0
      Prob.nState = double(k==0);
      Prob.Mode   = 2;
      r_k = nlp_r(x_k, Prob, varargin{:} );     % Function value
      if isempty(r_k)
         fprintf('The residual is empty. No NLLS problem!!!\n');
         Result.Iter=0;
         Result.IterGN=0;
         Result.ExitFlag=1;
         Result.ExitText='The residual is empty. No NLLS problem';
         Result.Inform=105;
         Result.x_0=x_k;
         Result.f_k=NaN;
         Result.x_k=x_k;
         Result=endSolve(Prob,Result);
         return
      end
      if any(isnan(r_k)) | any(isinf(r_k))
         fprintf('The residual contains NaN or Inf values at x_0\n');
         fprintf('Restart with a Prob.x_0 that is finite computable\n');
         Result.Iter=0;
         Result.IterGN=0;
         Result.ExitFlag=1;
         Result.ExitText='NaN or Inf values in the residual at start';
         Result.Inform=105;
         Result.x_0=x_k;
         Result.f_k=NaN;
         Result.x_k=x_k;
         Result=endSolve(Prob,Result);
         return
      end
      M = length(r_k);
      Prob.Mode   = 1;
      if alg ~= 7
         J_k = nlp_J(x_k, Prob, varargin{:} );     % Jacobian
      else
          if isfield(Prob.LS,'J_k')
             J_k = Prob.LS.J_k;
          else
             J_k = nlp_J(x_k, Prob, varargin{:} ); % Estimate the jacobian
                                                   % if not given.
          end
      end
      if any(any(isnan(J_k))) | any(any(isinf(J_k)))
         fprintf('The Jacobian contains NaN or Inf values at x_0\n');
         fprintf('Restart with a Prob.x_0 that is finite computable\n');
         Result.Iter=0;
         Result.IterGN=0;
         Result.ExitFlag=1;
         Result.ExitText='NaN or Inf values in the Jacobian at start';
         Result.Inform=105;
         Result.x_0=x_k;
         Result.f_k=NaN;
         Result.x_k=x_k;
         Result=endSolve(Prob,Result);
         return
      end
      f_k = 0.5 * (r_k' * r_k);
      g_k = J_k' * r_k;
      if k==0, gProj = g_k(idx); end
      if k > 0
         if m > 0
            Prob.Mode   = 2;
            c_k  = nlp_c( x_k, Prob, varargin{:} );
            Prob.Mode   = 1;
            dc_k = nlp_dc(x_k, Prob, varargin{:} );
         end
         if mA > 0, Ax=Prob.A*x_k; end
         if mTot > 0
            z=[Ax(:);c_k(:)];
            zErr = max(0,max(bl-z,z-bu));
            zSum = sum(zErr);
         end
      end
      if k==0
         Result.f_0=f_k;
         Result.x_0=x_k;
         JJ = J_k'*J_k;
         if alg==3
            H_k= JJ + A_k;
         end
         if alg == 5 | alg == 6
            % If JJ not full rank
            % A_k = 0.1*sqrt(f_k)*speye(Prob.N,Prob.N);
            % otherwise A_k = 0
            % B_k = JJ + A_k;
            B_k = JJ;
            H_k = JJ;
         end
         if PriLev > 4
            if ~isempty(Ud2r)
               H2 = nlp_d2r(x_k,Prob,r_k,J_k,varargin{:});
               if ~isempty(H2)
                  fprintf('Eigenvalues at start for Analytic Hessian\n');
                  xprinte(eig(JJ+H2),'eig:');
               end
            end
         end
      end
      Prob.nState = 0;
   end
   Result.x_k=x_k;
   Result.g_k=g_k;
   Result.c_k=c_k;
   Result.cJac=dc_k;

   %[KT, KTT, h_k, cErr] = KuhnTucker(Prob, Result,...
   %   bTol, zTol, xTol, x_k, x_L, x_U, g_k, z, bl, bu, xEqual, zEqual);

   if IterPrint
      fprintf(' %5d  %7d   %20.17f %11.7e %10.6f %17.7e\n', ...
              k, n_r, f_k, norm(p_full), alpha, norm(gProj));
   end

   itPrint(PriLev,x_k,f_k,g_k,r_k,J_k,k,pRank,pCond,nract,set,zSum,KT,KTT);

   release=0;
   releasecon=0;
   relcon = -1; % Constraint released in this iteration
   xR=[]; xRp=[];  % Var# released and gradient value


   if mi <= 0 % If no inequality constraints
      %
      % Check if any active variable shall be released
      %
      xR=[]; xRp=[];  % Var# released and gradient value
      % Negative grad on upper bound, positive on lower
      % ALT
      %lam_L = -set(:) .* g_k;
      % END ALT
      lam_L = -set(:) .* g_k;
      [minL ixL]= sort(lam_L(b_idx)./max(1,abs(x_k(b_idx))));

      % Dangerous to use 1st order estimate. May ==> infinite loop
      if nract > 0 & alpha0cnt < 3
         if PriLev > 2
            fprintf('1st order Lagrange multiplier estimate. ');
            xprinti(b_idx,'Var#:');
            xprinte(lam_L(b_idx),'Lam:');
         end
         % ALT
         %for i=1:nract
         %   j=b_idx(i);
         % END ALT

         for ii=1:nract
            % Take the variable with the most negative Lagrange multiplier 1st
            i=ixL(ii);
            j=b_idx(i);
            if nract < n
               % If lam_L > 0 optimal. Change sign in test
               % Release vars with negative Lagrange multiplier if they not
               % just hit a bound.
               % Also try a gradient step if releasing a vars result in 0 step.
                if (double((set_0(j)==set(j)) & (x_U(j)~=x_L(j)))*lam_L(j) < -xTol) |...
                   ((set_0(j)==0)  & ...
                   double((alpha==0) & (x_U(j)~=x_L(j)))*lam_L(j) < -xTol)
                   xR=[xR;j];
                   xRp=[xRp;-g_k(j)];
                   set(j)=0;
                   release=1;
                   % ALT, no break
                   break;
                end
             else % No freedom to move. Must check directly if optimum
                if (double(x_U(j)~=x_L(j)))*lam_L(j) < -xTol
                   set(j)=0;
                   release=1;
                   % ALT, no break
                   break;
                end
             end
         end
      end
      lambda = lam_L;
   else % There are inequality constraints

      if k==0 % Set up active constraint set W_kset
         W_km1set = W_kset;
         % Include all equalities in working set
         W_kset = [ones(me,1);zeros(mi,1)];
         actIcset = zeros(mc,1);
         con = b-A*x_k;
         if any(con > bTol)
            if PriLev > 0
               fprintf('\n ERROR ! x_k not feasible')
               fprintf(', constraints violated:\n')
               xprinti(find(con > bTol),'index');
               xprinte(x_k,' x_k ');
            end
            Result.x_k=x_k;

            Result.Iter=k;
            Result.IterGN=IterGN;
            Result.ExitFlag=1;
            Result.ExitText='x_k not feasible, constraint violated';
            Result.Inform=104;

            Result.f_k=f_k;
            Result.g_k=g_k;
            if ~isempty(H2)
               Result.H_k=JJ+H2;
            end
            Result.B_k=B_k;
            Result.v_k = LagMult(Prob,Result);
            %Result.v_k=[lam_L(:);lambda(:)];
            Result.r_k=r_k;
            Result.J_k=J_k;
            Result=endSolve(Prob,Result);
            return
         end
         ix=find(abs(con(me+1:mc)) <= bTol); % Active
         nractcon=min(length(ix),max(n,me)-(me+nract));
         if nractcon > 0
            W_kset(me+ix(1:nractcon))=1;
            actIcset(me+ix(1:nractcon))=1;
         end
      end

      ConTest(b, A, x_k, me, bTol, PriLev)

      W_k = find(W_kset); % Working set
      actIc = find(actIcset); % Set of active inequality constraints

      % If any changes in working set, update Aw_k
      %if any(W_kset ~= W_km1set)
      %   Aw_k = A(W_k,:);
      %end


      if PriLev > 2
         if me+nractcon > 0
            fprintf('Active constraints (including equalities) = %3.0f\n',...
                     nractcon+me);
            xprint(W_kset,'W_kset:',' %2.0f',25);
         end
      end

      if releaseall & ~isempty(xRR) & (alpha < bTol) % What tolerance to use ?
         % Releasing all vars with neg Lagr. mult. was a mistake
         for i = 1:length(xRR)
            if nract+nractcon < n
               set(xRR(i))=xRRlogic(i);
               nract = nract + 1;
               if PriLev > 2
                  fprintf('Reactivating variable %d\n',xRR(i));
               end
            end
         end
         nract = sum (set~=0);
         b_idx = find(set);
         idx = find(~set);
         releaseall = 0;
      else
         releaseall = 1;
      end

      if NOT_release_all % mbk debug, release max 1 var.
         releaseall=0;
      end

      % Check if any active constraint or variable is to be released
      xRR=[];xRRp=[];
      if (nractcon + nract > 0) % & pNorm < 1e-6
         % Estimate Lagrange multiplier for active constraints by solving
         % the overdetermined system: B'*lambda = g_k, where B is the active
         % inequality rows in A(W_k,:) augmented with simple bounds

         % Variables active on lower bound
         iL=find(set==-1);
         mL=length(iL);

         % Variables active on upper bound
         iU=find(set==1);
         mU=length(iU);

         % fixbnd is used together with lambda to determine if test
         % for convergence should be performed.

         fixbnd = [zeros(nractcon,1);xEqual(iL);xEqual(iU)];

         [Q, R, E, pRank] = ComputeQR(full([A(W_k,:);...
             sparse(1:mL,iL,ones(mL,1),mL,n); ...
             sparse(1:mU,iU,-ones(mU,1),mU,n)]'), epsRank);

         v_k = tomsol(6, Q, R, E, pRank, g_k);
         lambda = v_k(me+1:length(v_k));


         if PriLev > 2
            fprintf('Lagrange multiplier estimate: ');
            xprinte(v_k,'v_k:');
         end

         % Don't care about negative lambdas corresponding to fixed variables
         % when determining which constraint or which variable to release.
         lambda_hat = lambda;
         lambda_hat(find(fixbnd)) = 1;

         % Min cost rule, should Blands rule be used ?
         [lambdamin q_hat]=min(lambda_hat);

         %if lambdamin < -1E-8
         if lambdamin < -bTol

            if q_hat <= nractcon % Release constraint
               qq = actIc(q_hat);
               if qq ~= actcon; % Do not release just activated constraint
                  W_kset(qq) = 0;
                  actIcset(qq) = 0;
                  releasecon = 1;
                  relcon = qq;
                  % NOT USED ! Constraint relcon shall not be considered in
                               % computation of alpha_c if not violation
                  %relcon = -1;
               end
            elseif (q_hat > nractcon) & (q_hat <= nractcon+mL) & ~releaseall
               % Release var. from x_L
               jj = iL(q_hat-nractcon);
               if (jj ~= actvar) & (x_L(jj) ~= x_U(jj))
                  set( jj ) = 0;
                  xR = jj;
                  xRp = -g_k(jj);
                  %xR=[xR;jj];     % If releasing more than one variable
                  %xRp=[xRp;-g_k(jj)];
                  release = 1;
               end
            elseif (q_hat > nractcon+mL) & ~releaseall
               % Release variable from x_U
               jj = iU(q_hat-nractcon-mL);
               if (jj ~= actvar) & (x_L(jj) ~= x_U(jj))
                  set( jj ) = 0;
                  xR = jj;
                  xRp = -g_k(jj);
                  %xR=[xR;jj];     % If releasing more than one variable
                  %xRp=[xRp;-g_k(jj)];
                  release = 1;
               end
            end
            if releaseall % Release all vars with negative lagrange multiplier
               %xRR=[];xRRp=[]; This row is moved
               iLU = [iL;iU];
               lam = lambda(nractcon+1:length(lambda));
               for i = 1:length(iLU)
                  %if (iLU(i)~=actvar)&(lam(i)<-1E-8)&(x_L(iLU(i))~=x_U(iLU(i)))
                  if (iLU(i)~=actvar)&(lam(i)<-bTol)&(x_L(iLU(i))~=x_U(iLU(i)))
                     xRR  = [xRR ;iLU(i)];
                     xRRp = [xRRp;-g_k(xRR)];
                  end
               end
               if ~isempty(xRR)
                  xRRlogic = set(xRR);
                  set(xRR) = 0;
                  xR=xRR; xRp=xRRp;
                  release = 1;
               end
               %xRR = iLU(lam < -1E-8);
               %xRR = iLU(lam < -bTol);
               %set(xRR) = 0;
               %release = 1;
               %xRRlogic = set(xRR);
               %%iLU
               %%lam
               %%xRR
            end
         end
      else
         fixbnd = 0; % Used together with lambda to determine if test
                     % for convergence should be performed.
      end
      actcon = -1;
      actvar = -1;
   end

%
% Check convergence conditions, if not releasing variables or constraints
%
   Inform = 0;
   if release | releasecon
      if release
         nract = sum (set~=0);
         b_idx = find(set);
         % Index vector for non active variables, the free vars
         idx = find(~set);
         if PriLev > 2
            fprintf('Release variable from bound. ');
            xprinti(xR,'var#:');
         end
      end
      if releasecon
         nractcon = sum(actIcset~=0);
         W_k = find(W_kset);
         actIc = find(actIcset);
         if PriLev > 2
            fprintf('Constraint %d is released\n',qq);
         end
      end
   end

   %
   % Update of Z, Null space basis
   %

   if LCflag & ( any(set_0 ~= set) | any(W_kset ~= W_km1set) | k==0 )
      if isempty(W_k) | isempty(idx)
         Z=speye(length(idx));
         Arank=0;
      else
         % Permutation of rows (constraints) in QR, to get null space last in Q
         [Q, R, E, Arank] = ComputeQR(A(W_k,idx)', epsRank);
         Z=sparse(Q(:,Arank+1:length(idx)));
      end
   end
   if mi > 0
      W_km1set=W_kset;
      set_0 = set;
   end

   if ~(release | releasecon)% | fRed0cnt > LowIts % Check convergence criteria
      %if all(lam_L >=-xTol) % Only check if possible optimum

      if mi <= 0
         fixbnd = xEqual;
      else

      end
      % What tolerance to use ????
      % Only check if possible optimum

      %if all(lambda >=-1E-8 | fixbnd)  % | fRed0cnt > LowIts
      if all(lambda >=-bTol | fixbnd)
         if (me+mi) > 0
            gProj = full(Z'*g_k(idx));
         else
            gProj = g_k(idx);
         end

         % Convergence criteria 1
         if ~isempty(gProj)
            gOK=Robust*max(abs(gProj) * max(max(abs(x_k(idx)),size_x)))...
                  <= 10000*optParam.eps_g * max(abs(f_k),size_f);
         else
            gOK=1;
         end
         %if max(abs(x_k-x_km1)./max(abs(x_k),size_x)) <= eps_x   & ...
         if max(abs(x_k-x_km1)./max(1,abs(x_k))) <= eps_x   & ...
            (alphaMax > 1E-14) & gOK & TESTONX
            %& gOK & TESTONX

            if PriLev >= 1
               disp('*** Convergence 1, Iteration points are close ***');
            end
            Inform=Inform+1;
            flag=0;
            stop = 1;
         end

         % Check for convergence on projected gradient
         % The dimension of x_k(idx) is not OK
         % optParam.eps_g;    % Gradient convergence tolerance
         %if max(abs(gProj).* max(abs(x_k(idx)),size_x)) <= ...

         if ~isempty(idx) & alg~= 7
            if max(abs(gProj) * max(max(abs(x_k(idx)),size_x))) <= ...
                   optParam.eps_g * max(abs(f_k),size_f) & ...
               Robust*max(abs(x_k-x_km1)./max(abs(x_k),size_x)) <= 100*eps_x
               if PriLev >= 1
                  disp('*** Convergence 2, Projected gradient small ***');
               end
               Inform=Inform+2;
               flag=0;
               stop = 1;
            end
         end

         % Simple gradient test, for Wang-Li-Qi special test
         %if norm(gProj) <= optParam.eps_g
         %   if PriLev >= 1
         %      disp('*** Convergence 2, Projected gradient small ***');
         %   end
         %   Inform=Inform+2;
         %   flag=0;
         %   stop = 1;
         %end

         if abs(f_k-fLow) <= eps_absf*max(1,fLow) % Convergence criteria 3
            if PriLev >= 1
               disp('*** Convergence 3, Function value close to Prob.f_Low');
            end
            flag=0;
            Inform = Inform + 4;
            stop = 1;
         end

         if fRed0cnt > LowIts 	      % Convergence criteria 4
            if PriLev >= 1
               fprintf('*** Convergence 4, Relative function value reduction')
               fprintf(' low for %d iterations ***\n',LowIts);
            end
            flag=0;
            Inform = Inform + 8;
            stop = 1;
         end
         %HKH
         %if 0
         %if f_red < 1E-14*max(1,f_k)   % Convergence criteria 5
         %if f_red < 1E-30*max(1,f_k) & alg~=7   % Convergence criteria 5
         if f_red < 1E-30*max(1,f_k) & alg==7   % Convergence criteria 5
            if Inform == 0 | NewJ == 0
               J_k = nlp_J(x_k, Prob, varargin{:} );     % Jacobian
               B_k = J_k;
               NewJ = 1;
               %HKH'recompute J_k'
               %Inform
               Inform=0;
               stop=0;
            else
               flag=0;
               Inform = Inform + 16;
               stop = 1;
            end
         elseif f_red < 1E-30*max(1,f_k)
         %if f_red < 1E-30*max(1,f_k)
            if PriLev >= 1
               disp('*** Convergence 5, Small change in f_k ***');
            end
            if sigma > 0.1
               %sigma = sigma*0.1;
               sigma = 0.1;
               Prob.LineParam.sigma = sigma;
            elseif sigma > 0.01
               sigma = 0.01;
               Prob.LineParam.sigma = sigma;
            %elseif alg ~= 7
            else
               flag=0;
               Inform = Inform + 16;
               stop = 1;
            end
         end
         %end
      end

%
% CUT
%

      if nract==n
         stop=1;
         flag=0;
         Inform=Inform+32;
         if PriLev >= 1
            disp('*** Convergence 5, Local min with all variables on bounds!');
            disp('*** No check has been done on the gradient.');
         end
      end
   end

   if k >= MaxIter		% Stop criteria 1
      if PriLev >= 1
         disp('*** STOP 1! Max no of iterations reached ***');
      end
      Result.ExitText='Maximal number of iterations reached';
      if Inform >0 & Inform < 100, flag=0;
      else
          flag=1; Inform=101;
      end
      stop = 1;
   end

   if f_k < fLow     % Stop criteria 2
      if PriLev >= 1
         disp('*** STOP 2, Function value below given estimate ***');
         disp('*** Restart with lower Prob.f_Low. ***');
         disp('*** if minimum not reached ***');
      end
      Result.ExitText=str2mat('Function value below Prob.f_Low. ' ...
             ,'Set lower if not minimum');
      if Inform >0 & Inform < 100, flag=0;
      else
          flag=0; Inform=102;
      end
      stop = 1;
   end
   if pRank == 0 & k == 1     % Convergence criteria 99
      if PriLev >= 1
         fprintf('*** Convergence 99, Residual independent of x, Jacobian == 0')
         fprintf('\n');
      end
      flag=0;
      Inform = 99;
      stop = 1;
   end

   if stop			% Show results before return
      if PriLev >= 5   % Results are shown in lsrun instead.
         fprintf('==================================================\n');
         fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);
         if ~(isinf(KT) | isnan(KT))
            fprintf('Kuhn-Tucker error                  %30.20f\n',max(KT));
         end

         xprint(x_k,'x_k:');

         xprint(set,'Set:',' %2.0f',25);

         xprinte(g_k,'g_k:');

         xprint(r_k,'r_k:');

         if ~isempty(H2)
            fprintf('The Analytic Hessian H_k\n');
            H2 = nlp_d2r(x_k,Prob,r_k,J_k,varargin{:});
            mPrint(JJ+H2,'H_k',' %13.6e',5);
         else
            fprintf('The Approximative Quasi-Newton Hessian B_k\n');
            mPrint(B_k,'B_k',' %13.6e',5);
         end


         if optParam.wait, pause; end;

      end
      Result.Iter=k;
      Result.IterGN=IterGN;
      Result.ExitFlag=flag;
      Result.Inform=Inform;
      switch  Inform
         case 1
           Text = 'Iteration points are close';
         case 2
           Text = 'Projected gradient small';
         case 3
           Text = 'Iterations close, small projected gradient';
         case 4
           Text = 'Function value close to 0';
         case 5
           Text = 'Iteration close, f(x) close to 0';
         case 6
           Text = 'Projected gradient small, f(x) close to 0';
         case 7
           Text = str2mat('Iterations close, small projected gradient' ...
                  ,'f(x) close to 0');
         case 8
           Text = 'Relative f(x) reduction low for LowIts iterations';
         case 11
           Text = 'Rel f(x) red. low for LowIts iter. Close Iters';
         case 16
           Text = 'Small Relative f(x) reduction';
         case 17
           Text = 'Close iteration points, Small relative f(x) reduction';
         case 18
           Text = 'Small gradient, Small relative f(x) reduction';
         case 32
           Text = 'Local minimum with all variables on bounds';
         case 99
           Text = 'The residual is independent of x. The Jacobian is 0';
         otherwise
           Text = ['Unknown Inform value ' num2str(Inform)];

      end
      if Inform < 100
         Result.ExitText=str2mat('Optimal solution found',Text);
      end

      Result.f_k=f_k;
      Result.g_k=g_k;
      Result.x_k=x_k;
      Result.v_k = LagMult(Prob,Result);
   %  Result.v_k=[lam_L(:);lambda(:)];
      Result.r_k=r_k;
      Result.J_k=J_k;
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=double(x_k==x_L)+2*double(x_k==x_U);
      if ~isempty(Prob.A)
         Ax=Prob.A*x_k;
         Result.bState=double(Ax<=b_L)+2*double(Ax>=b_U);
      else
         Result.bState=[];
      end
      Result=endSolve(Prob,Result);
      return;
   end
%
% Compute the Hessian
%
   if PriLev > 4
      % Print the Hessian, if possible
      if ~isempty(Ud2r)
         H2 = nlp_d2r(x_k,Prob,r_k,J_k,varargin{:});
         if ~isempty(H2)
            fprintf('The Analytic Hessian H_k\n');
            mPrint(JJ+H2,'H_k',' %13.6e',5);
         end
      end
   end


% Print the Quasi-Newton approximation to the Hessian

   if PriLev > 5 & ~isempty(B_k)
      fprintf('The Approximative Quasi-Newton Hessian B_k\n');
      mPrint(B_k,'B_k',' %13.6e',5);
   end


%
% Determine search direction p
%
   release=1;
   while release & ( ~(alpha==0 & ~isempty(xR)) | LCflag )

      [p, pRank, pCond]=SearchStep(r_k, J_k, g_k, B_k, Z, GN, LCflag,...
                        idx, epsRank, Method, n, nract, PriLev, LargeScale);
      release=0;
   end
   IterGN = IterGN + GN;

   p_full=zeros(n,1);
   if alpha==0 & ~isempty(xR) & ~LCflag
      % Search in the coordinate direction of the variables to be released.
      % Last full step was a total failure
      if LCflag
         % Projected gradient
         p_full(xR)=full(Z'*xRp);
      else
         p_full(xR)=xRp;
      end
      p=p_full(idx);
   elseif LCflag & isempty(Z)
      % Not sure of this case, cls_prob #34, p = 0, pRank = 0 set
   elseif pRank == 0
      % Residual is independent of x, J == 0
      p=zeros(length(idx),1);
   elseif length(p) < length(idx)
      % Search step in subspace computed. Put into full space vector p_full
      p = full(Z*p);
      p_full(idx) = p;
   else
      p_full(idx) = p;
   end
   pNorm=norm(p);

%
% Find length of search step, alpha
%                                      1
   TESTONX=1;
   fp_0=g_k(idx)'*p;
   if k == 0
      alpha_1 = LineParam.InitStepLength;
   elseif pNorm==0
      alpha_1 = LineParam.InitStepLength;
   else
      df = max(f_km1-f_k, 10*eps_x);
      if fp_0 == 0
         if PriLev > 0
            disp('gn: ZERO fp_0!!');
            xprinte(x_k,'x_k');
            xprinte(g_k,'g_k');
            xprinte(g_k(idx),'g_k(idx)');
            xprinte(p,'p');
         end
      end
      alpha_1 = min(LineParam.InitStepLength, -2*df/fp_0);
      if alpha_1 < 0
         if PriLev > 0
            fprintf('ALARM: No descent direction in iter %d. ',k)
         end
         if fp_0 < -1E-10
            alpha_1 = -alpha_1;
            p_full=-p_full;
            p=-p;
            if PriLev > 0
               fprintf('Change sign on p and try\n')
               fprintf('Function reduction: %15.7e  ',df);
               fprintf('Directed derivative: %15.7e\n',fp_0);
            end
         else
            TESTONX=0;
            if isempty(Z)
               p=-g_k(idx);
            else
               p=full(-Z*Z'*g_k(idx));
            end
            p_full(idx)=p;
            fp_0=g_k(idx)'*p;
            alpha_1 = min(LineParam.InitStepLength, -2*df/fp_0);
            if PriLev > 0
               fprintf('Use negative gradient step\n')
               fprintf('Function reduction: %15.7e  ',-2*df/fp_0);
               fprintf('Directed derivative: %15.7e\n',fp_0);
            end
         end
      end
   end	% (if)

   if PriLev > 2
      fprintf('Directed derivative: %15.7e\n',full(fp_0));
      %fprintf('Best alpha step estimated: max(0.001,%15.7e)\n',alpha_1);
      fprintf('Best alpha step estimated: max(0.5,%15.7e)\n',alpha_1);
   end


   % HKH OBZ. Should we test on pNorm >= eps_x instead

   if pNorm > 0   % Only line search if non zero step
      %alpha_1=max(1E-3,alpha_1);   % Safe guard to avoid too small step
      alpha_1=max(0.5,alpha_1);   % Safe guard to avoid too small step

      %
      % Line search
      %
      if PriLev > 2
         xprinte(p_full,'p:  ');
         if PriLev > 3 | sum(set) > 0
            fprintf('Active set:\n');
            xprint(set,'Set:',' %2.0f',25);
         end
         if PriLev > 3
            if alpha_1 < 0.99
               fprintf('Line search start %10.6f\n',alpha_1)
            end
         end
      end

      Step = Inf*ones(length(idx),1);
      L=find(~isinf(x_L(idx)) & p < -xTol * max(1,abs(x_k(idx))));
      U=find(~isinf(x_U(idx)) & p >  xTol * max(1,abs(x_k(idx))));
      if ~isempty(L)
         Step(L)=(x_L(idx(L))-x_k(idx(L)))./p(L);
      end
      if ~isempty(U)
         Step(U)=(x_U(idx(U))-x_k(idx(U)))./p(U);
      end
      [alphaMax xvarB]=min(Step);
      if isempty(L)
         xvartype=1;
      elseif any(xvarB==L)
         xvartype=-1;
      else
         xvartype=1;
      end
      xvar=idx(xvarB);
      alpha0min=min(1E20,alphaMax);

      pScale=0;

      if alphaMax <= 0
         if PriLev > 1
            fprintf('No move possible!! Local minimum');
            fprintf(' - Var %d limits. alphaMax %f.',xvar,alphaMax);
            fprintf('\n')
            fprintf(' Low %15.10e Upp %15.10e x %15.10e p %15.10e\n',...
                      x_L(xvar),x_U(xvar),x_k(xvar),p(xvarB));
         end
      elseif isinf(alphaMax)
         alphaMax=10;
      elseif alg==7 & ...
         ((alphaMax < 1 & pNorm > 1E3) | (alphaMax > 1E5 & pNorm < 1E-5))
         % Scale the step if Gauss Newton
         pScale=alphaMax/10;
         p=pScale*p;
         p_full=pScale*p_full;
         alphaMax=10;
         alpha_1 =LineParam.InitStepLength;
      end
      f_km1 = f_k;
      r_km1 = r_k;
      J_km1 = J_k;  % Save Jacobian
      g_km1 = g_k;  % Used by Fletcher-Xu

%*******************************************************************************
      if mi > 0 % If any inequality constraints
         alpha_c = 1E20;
         alpha_cmin = 1E20;
         % Determine max steplength using  nonactive inequality constraints
         NactIc = find(~actIcset); % Set of nonactive inequality constraints
         for i = 1:length(NactIc)
            if (NactIc(i) > me ) & ( A(NactIc(i),:)*p_full < 0 )
            %& ( NactIc(i) ~= relcon )
               alpha_c = ( b(NactIc(i))-A(NactIc(i),:)*x_k ) / ...
                         ( A(NactIc(i),:)*p_full );
               if alpha_c < alphaMax
                  alphaMax = alpha_c;
                  alpha_cmin = alpha_c;
                  rr = NactIc(i);
                  if PriLev > 2
                     fprintf('Constraint %d restricts alphaMax',NactIc(i));
                     fprintf(' in iteration %d',k);
                     fprintf(', alpha_c=%8.6f\n',alpha_c);
                  end
               end
            end
         end
         % relcon = -1; This row is moved
      end
%*******************************************************************************

      % OBZ HKH quick change for Huschens method. Limit the line search step
      if alg==3, alphaMax=min(alphaMax,1); end

      % Should we really restrict the step this hard?
      % if alg < 3,  alphaMax=min(alphaMax,2); end

      if alphaMax <= 1E-14   % % Too small value. Can not take any step
         alpha=0;
         if PriLev > 2
            fprintf('alphaMax too small %15.7e. No line search possible\n',...
            alphaMax);
         end
         f_red=0;
      else

         LineParam.InitStepLength = alpha_1;
         LineParam.r_k = r_k;
         LineParam.J_k = J_k;
         if LineAlg <= 1
            LineResult = LineSearch( 'nlp_r','nlp_J', x_k, p_full, f_k, g_k,...
                         LineParam,alphaMax,2, PriLev-3,Prob,varargin{:});
         else
            LineResult=Armijo(x_k,p_full,f_k,g_k,agFac,sigma,alphaMax,...
                              Prob,varargin{:});
         end


%if 0
%         alpha=LineResult.alphaVec;
%         disp('Plot LineSearch result')
%
%         if 1 | plotLine
             Prob.Mode   = 0;
%            LinePlot('nlp_r', x_k, p_full, f_k, g_k, LineParam, alphaMax,...
%                     2, alpha, Prob, varargin{:});
%         end
%end
%if 0
%
%         LineResult = eloLine( 'nlp_r','nlp_J', x_k, p_full, f_k, g_k,...
%         r_km1,J_km1,LineParam,alphaMax, 2, PriLev-3,Prob,varargin{:});
%
%
%         alpha=LineResult.alphaVec;
%
%         disp('Plot eloLine')
%         if 1 | plotLine
             Prob.Mode   = 0;
%            LinePlot('nlp_r', x_k, p_full, f_k, g_k, LineParam, alphaMax,...
%                     2, alpha, Prob, varargin{:});
%         end
%end
         alpha=LineResult.alpha;
         %if alpha==0 & fp_0 < -1E-5 & pNorm > 1E-5
         if alpha==0
            % Scales are making trouble
            alpha_z=1E-8;
            alphaVec=LineResult.alphaVec;
            r_k=LineResult.r_k;
            J_k=LineResult.J_k;

            if alg==4
               alphaTest=alphaMax*1E-9*0.1.^[1:7];
            else
               if alphaMax < 1E-8
                  alphaTest=[alphaMax,alphaMax*[0.1.^[1:16]]];
               else
                  alphaTest=alphaMax*[1E-7*0.1.^[1:12]];
                  %alphaTest=alphaMax*[1E-9*0.1.^[1:12],0.1.^[9:-1:1]];
               end
            end
            alpha_z=1E-8;
            alpha_z=alphaMax;
            fLast=Inf;
            for ii=1:length(alphaTest);
                alpha_z=alphaTest(ii);
                Prob.Mode   = 0;
                r_z = nlp_r(x_k+alpha_z*p_full, Prob, varargin{:} );
                f_z = 0.5 * (r_z'*r_z);
                alphaVec=[alphaVec,alpha_z];
                if f_z < f_k
                   f_k=f_z;
                   r_k=r_z;
                   alpha=alpha_z;
                elseif f_z >= fLast & alpha > 0
                   break;
                end
                fLast=f_z;
            end
            if alpha > 0
               TESTONX=0;
               Prob.Mode   = 1;
               J_k = nlp_J(x_k, Prob, varargin{:} );     % Jacobian
               g_k = J_k' * r_k;
            end
         else
            f_k=LineResult.f_alpha;
            g_k=LineResult.g_alpha;
            alphaVec=LineResult.alphaVec;
            r_k=LineResult.r_k;
            if ~isempty(LineResult.J_k) & alpha > 1E-14
               J_k=LineResult.J_k;
            end
         end
         f_red=f_km1-f_k;
         if PriLev > 2
            fprintf('Line search: alphaMax=%10.4e.',alphaMax);
            fprintf(' Func.red=%10.4e. ',f_red);
            if f_k > 0
               fprintf(' Rel.red=%10.4e. ',f_red/f_k);
            end
            if alpha > 1E-5
               fprintf(' Step length=%12.6f. ITER %d\n',alpha,length(alphaVec));
            else
               fprintf(' Step length=%12.6e. ITER %d\n',alpha,length(alphaVec));
            end
         end
      end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      % Update
      %
      x_km1 = x_k;
      x_k = x_k + alpha*p_full;

      if alpha > 1E-30 | ~TESTONX
         JJ = J_k'*J_k;
         if m > 0
            Prob.Mode = 2;
            c_k  = nlp_c( x_k, Prob, varargin{:});
            Prob.Mode = 1;
            dc_k = nlp_dc(x_k, Prob, varargin{:});
         end
         if mA > 0, Ax=Prob.A*x_k; end
         if mTot > 0
            z=[Ax(:);c_k(:)];
            zErr = max(0,max(bl-z,z-bu));
            zSum = sum(zErr);
         end
      end

      % Determine if any constraint or variable is to be active

      if mi > 0
         if abs(alpha - alpha0min) < xTol % Variable shall be activated
            %set_0 = set;
            set(xvar) = xvartype;
            nract = nract +1;
            actvar = xvar;
            if PriLev > 2
               fprintf('Variable %d is activated\n',xvar);
            end
         elseif abs(alpha - alpha_cmin) < bTol % Constraint shall be activated
            %W_km1set = W_kset;
            W_kset(rr)=1;
            actIcset(rr) = 1;
            actcon = rr;
            nractcon = nractcon +1;
            if PriLev > 2
               fprintf('Constraint %d is activated\n',rr);
            end
         end
      end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if alpha < 1E-30
         alpha0cnt=alpha0cnt+1;
      else
         alpha0cnt=0;
      end
      if f_red <= 100 * eps * max(f_k,size_f)
         fRed0cnt=fRed0cnt+1;
      else
         fRed0cnt=0;
      end

      if optParam.wait, pause; end;

      if CURV~=0, xLast=x_km1; end

      if alg==4
         B_k=JJ;
      else
         % Check if we moved from x_k
         if alpha > 100 * eps & pNorm > 0
            [B_k, GN, A_k] = SafeUpdate( alg, alpha, p_full, eps_x, f_k, ...
                         f_km1, J_k, J_km1, r_k, r_km1, B_k, JJ, g_k, g_km1, ...
                         A_k, PriLev, k, pCond, pRank, GN);
            if alg==7
               J_k = B_k;
               g_k = J_k'*r_k;
            end
         elseif alg==7
            J_k = nlp_J(x_k,Prob);
            g_k = J_k'*r_k;
         end

      end % Check if we moved from x_k

   else % Only line search if non zero step
      x_km1=x_k;
      if optParam.wait, pause; end
   end % Only line search if non zero step

   % OBZ! Should we count a zero step as an iteration?
   %      Or just a change of active set?
   k = k+1;
end	% (while)

% =========================================
function [B_k, GN, A_k] = SafeUpdate(alg, alpha, p_full, eps_x, f_k, f_km1,...
   J_k, J_km1, r_k, r_km1, B_k, JJ, g_k, g_km1, A_k, PriLev, k, pCond,pRank, GN)
% =========================================

if alg <= 1 % Fletcher-Xu
   z=alpha*p_full; % Check if we have converged (or short alpha step)
   if ((f_km1-f_k) >= 0.2*f_km1 | norm(z) <= eps_x)
      GN=1;
      B_k=JJ;
   else            % BFGS update. SAFE GUARDED!
      GN=0;
      Bz=B_k*z;
      zBz=z'*Bz;

      y=JJ*z+(J_k-J_km1)'*r_k;
      y_old=g_k-g_km1;
      zy=z'*y;
      zy_old=z'*y_old;
      if zy < 0.01 * zy_old
         w=y_old;
         zw=zy_old;
      else
         w=y;
         zw=zy;
      end
      if PriLev >= 3
         fprintf('BFGS update iter %3.0f\n',k);
         fprintf('zw %10.6e zBz %10.6e \n',zw,zBz)
      end
      if zw < 1000 * eps | zBz < 1000 * eps
         if PriLev > 2
            %fprintf('BFGS step dangerous! Use Gauss-Newton\n',k);
            fprintf('BFGS step dangerous! Restart with Identity\n');
         end
         if PriLev > 3
            fprintf('max(w)/zw %20.10e\n',max(w)/zw);
            fprintf('min(w)/zw %20.10e\n',min(w)/zw);
            fprintf('max(Bz)/zBz %20.10e\n',max(Bz)/zBz);
            fprintf('min(Bz)/zBz %20.10e\n',min(Bz)/zBz);
         end

         %GN=1;
         %B_k=JJ;
         if PriLev > 1
            disp('BFGS RESTART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
         end
         GN=0;
         B_k=speye(size(JJ));
      else
         B_k=B_k+(w/zw)*w'-(Bz/zBz)*Bz';
      end
      %if pCond > 1E7 & GN==0
      %   GN=1;
      %end
      if PriLev > 3
         fprintf('B_k matrix after BFGS update\n');
         PrintMatrix(B_k,'B_k:');
         fprintf('B_k matrix eigenvalues:\n');
         xprinte(eig(B_k),'eig:');
      end
   end
end

if alg == 2 % HKH:s safeguard to  Al-Baali-Fletcher hybrid method
   z=alpha*p_full; % Check if we have converged (or short alpha step)
   if (f_km1-f_k) >= 0.2*f_km1 | norm(z) <= eps_x
      GN=1;
      B_k=JJ;
   else % BFGS update. SAFE GUARDED!
      GN=0;
      y=JJ*z+(J_k-J_km1)'*r_k;
      Bz=B_k*z;
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
         fprintf('BFGS update iter %3.0f\n',k);
         fprintf('zw %10.6e zBz %10.6e \n',zw,zBz)
      end
      if zw < 1000 * eps | zBz < 1000 * eps
         if PriLev > 2
            fprintf('BFGS step dangerous! Use Gauss-Newton\n');
         end
         %GN=1;
         %B_k=JJ;
         if PriLev > 1
            disp('BFGS RESTART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
         end
         GN=0;
         B_k=speye(size(JJ));
      else
         B_k=sparse(B_k+(w/zw)*w'-(Bz/zBz)*Bz');
      end
      if PriLev > 3
         fprintf('B_k matrix after BFGS update\n');
         PrintMatrix(B_k,'B_k:');
         fprintf('B_k matrix eigenvalues:\n');
         xprinte(eig(B_k),'eig:');
      end
   end
end
if alg == 3 % Huschens TSSM algorithm. Always BFGS type of update.
   GN=0;   % From now on take QuasiNewton steps
   %z=alpha*p_full;             % The step. TSSM: var. s
   z=p_full;                    % The step. TSSM: var. s
   % HERE MAYBE SAFE GUARDING SHOULD BE INCLUDED, test on f_km1*2
   y_cr=(J_k-J_km1)'*r_k/sqrt(2*f_km1);    % TSSM: y^#
   norm_r_k=sqrt(2*f_k);                   % TSSM: |R(x+)|
   y=JJ*z+norm_r_k*y_cr;                   % TSSM: y
   Bs=JJ+norm_r_k*A_k;                     % TSSM: B^s
   tau=1;
   Bs_z=Bs*z;
   q=z'*Bs_z;
   yz=y'*z;  % This scalar product may become slightly < 0
   if q > 0 & yz > 0    % Safe guarding close to optimum
      % sig may become imaginary. Use real() or test yz > 0
      sig=sqrt(yz/q);
   else
      if PriLev > 1
         fprintf(' Huschens: Quote <= zero:%30.20f\n',q);
         fprintf('Step length alpha %20.10f. ',alpha);
      end
      sig=0;
   end
   A_k = A_k+secUpdat(z,y_cr,A_k,y+tau*sig*Bs_z);% TSSM: A+
   B_k = JJ+norm_r_k*A_k;                        % TSSM: B+

   %if 0   % Check the theory on page 112 in SIAM -94 paper
   %   Bp=Bs+secUpdat(z,y,Bs,y+tau*sig*Bs_z);
   %   fprintf('Matrices should be equal:');
   %   %sum(Bp~=B_k)
   %   PrintMatrix(Bp-B_k,'Bp-B_k:') % Pretty equal, about 1E-10 error
   %   pause
   %end

   if PriLev > 2
      fprintf('B_k matrix after Huschens TSSM BFGS update.');
      fprintf(' Cond(B_k) %16.8e\n',cond(B_k));
   end
   if PriLev > 3
      PrintMatrix(B_k,'B_k:');
      fprintf('B_k matrix eigenvalues:\n')
      xprinte(eig(B_k),'eig:');
   end

end
if alg == 5 | alg == 6 % alg = 5: Wang, Li, Qi Structured MBFGS method
                       % alg = 6: Li-Fukushima MBFGS method
   z     = alpha*p_full; % Check if we have converged (or short alpha step)
   n     = length(z);
   gNorm = norm(g_k);
   JJz   = JJ*z;
   yC    = (J_k-J_km1)'*r_k;
   yHat  = JJz + yC;   %yHat = JJz + (J_k-J_km1)'*r_k;
   zz    = z'*z;
   yHatz = yHat'*z;
   if gNorm > 1
      beta = gNorm^0.01;
   else
      beta = gNorm*gNorm;
   end
   if yHatz >= 0 |  zz == 0
      tk = 1E-6*beta;
   else
      tk = beta-yHatz/zz;
   end
   if alg == 5
      % alg=5: y is y = JJ*z+(J_k'-J_km1')*r_k + tk*z;
      y  = yHat + tk*z;
   else
      % alg=6: y is y = (J_k'-J_km1')*r_k + tk*z;
      y  = yC   + tk*z;
   end
   % yC  = yC   + tk*z;   % yC = y_k^#
   Bz  = B_k*z;
   zBz = z'*Bz;
   zy  = z'*y;
   % v   = y + sqrt(zy/zBz)*Bz;

   if (f_km1-f_k) >= 0.2*f_km1 | norm(z) <= eps_x
      % If JJ is well-conditioned: (use pRank)
      % B_k = JJ;
      % else, add term on diagonal, factor * norm(r_k) * I
      % B_k = JJ + 0.1*sqrt(f_k)*speye(n,n);
      % Note - use sqrt(f_k) = sqrt(0.5*r_k'r_k), not e.g. |r_k|
      if pRank < n & GN == 1
         B_k = JJ + 0.1*sqrt(f_k)*speye(n,n);
      else
         B_k = JJ;
      end
      GN  = 1;
   else % MBFGS update.
      GN  = 0;
      % A_ks = B_k - JJ;
      % The following two secant updates should be the same
      % Up1 = secUpdat(z,y,B_k,v);
      % Up2 = secUpdat(z,yC,A_ks,v);
      % 'secant update check'
      % sum(Up1-Up2)

      % Check update. A_k used before new update
      % B3 = JJ + A_k + secUpdat(z,y,JJ+A_k,v); % Secant update

      % A_k = A_ks + secUpdat(z,yC,A_ks,v); % Secant update

      % Update using secUpdat formula
      % B_k = JJ + A_k;

      % Simple update formula
      B_k = B_k - Bz*Bz'/zBz + (y/zy)*y';

      if PriLev > 3
         fprintf('B_k matrix after BFGS update\n');
         PrintMatrix(B_k,'B_k:');
         fprintf('B_k matrix eigenvalues:\n');
         xprinte(eig(B_k),'eig:');
      end
   end
end
if alg==7 % Broydens method
   GN = 1;
   z  = alpha*p_full; % Check if we have converged (or short alpha step)
   s  = r_k-r_km1;
   zz = z'*z;
   if zz > 1000 * eps | norm(z) <= eps_x
      B_k = J_k +((s-J_k*z)/zz)*z';
   else
      B_k = J_k;
   end
end

% =========================================
function [alg, SolvAlg, B_k, A_k]=NameAlg(Prob)
% =========================================

alg=max(0,Prob.Solver.Alg);
if isempty(alg), alg=0; end

if alg==1
   SolvAlg='Fletcher-Xu Hybrid Method. Modified Gauss-Newton - BFGS';
elseif alg==2
   SolvAlg='Al-Baali - Fletcher Hybrid Method. Modified GN - BFGS';
elseif alg==3
   SolvAlg='Huschens TSSM algorithm';
elseif alg==4
   SolvAlg='Gauss Newton with Subspace Minimization';
elseif alg==5
   SolvAlg='Wang, Li, Qi Structured MBFGS method';
elseif alg==6
   SolvAlg='Li-Fukushima MBFGS method';
elseif alg==7
   SolvAlg='Broydens method';
else
   % Default algorithm
   alg=1;
   SolvAlg='Fletcher-Xu Hybrid Method. Modified Gauss-Newton - BFGS';
end
if Prob.LargeScale
   SolvAlg=[SolvAlg '. Iterative LSQR.'];
end

% Find initial Quasi-Newton matrix
n = Prob.N;
if alg<=3 | alg >= 5
   if ~isempty(Prob.optParam.QN_InitMatrix)
      B_k = Prob.optParam.QN_InitMatrix;   % B = User given matrix
      if size(B_k,1)~=n | size(B_k,2)~=n
         disp('Illegal dimensions on initial Quasi Newton matrix');
         disp('clsSolve is using identity matrix instead');
         B_k = speye(n);                % B = Use identity matrix
      end
   else
      B_k = speye(n); % B_k = Start Quasi-Newton with identity matrix
   end
else
   B_k = [];
end
if alg==3
   A_k=B_k;  % Initial value of Huschens update matrix
else
   A_k=[];
end


% =========================================
function [p, pRank, pCond]=SearchStep(r_k, J_k, g_k, B_k, Z, GN, LCflag,...
             idx, epsRank, Method, n, nract, PriLev, LargeScale)
% =========================================
pCond=[];
pRank=[];
if LCflag & isempty(Z)
   %fprintf('\n Z is empty\n');
   pRank=0;
   p = zeros(length(idx),1);
   return
end

if Method == 0 & GN  % Gauss-Newton with QR-Decomposition, and pivoting

   if LargeScale
      if isempty(Z)
         %p1 = qls(sparse(J_k(:,idx)), -r_k);
         %[p, pRank]  = ssqls(sparse(J_k(:,idx)), -r_k);
         [p, pRank]  = sls(sparse(J_k(:,idx)), -r_k);

%[ p, iStop, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
%      Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
%             WantStdErr, nOut, x_0 )


      else
         %p1 = qls(sparse(J_k(:,idx)*Z), -r_k);
         %[p, pRank]  = ssqls(sparse(J_k(:,idx)*Z), -r_k);
         [p, pRank]  = sls(sparse(J_k(:,idx)), -r_k);

      end
   else
      if isempty(Z)
         [Q, R, E, pRank] = ComputeQR(J_k(:,idx), epsRank);
      else
         [Q, R, E, pRank] = ComputeQR(J_k(:,idx)*Z, epsRank);
      end


      if pRank==0
         pCond=Inf;
      else
         pCond=abs(R(1,1))/max(1E-200,abs(R(pRank,pRank)));
      end
      if pRank == 0
         if any(isnan(g_k)) | any(any(isnan(B_k)))
            fprintf('ERROR!!! THE USER HAS RETURNED NAN ELEMENTS ')
            fprintf('IN THE RESIDUAL OR JACOBIAN\n')
         end
      end

      if isempty(Z)
         p = tomsol(6, Q, R, E, pRank, -r_k);
      else
         p = full(Z * tomsol(6, Q, R, E, pRank, -r_k));
      end
   end

elseif Method == 0 & ~GN  % BFGS, using QR-Decomposition, and pivoting

   if LargeScale
      if isempty(Z)
         %p1 = qls(sparse(B_k(idx,idx)), -g_k(idx));
         %[p, pRank] = ssqls(sparse(B_k(idx,idx)), -g_k(idx));
         [p, pRank]  = sls(sparse(B_k(idx,idx)), -g_k(idx));
      else
         %p1 = qls(sparse(Z'*B_k(idx,idx)*Z), -Z'*g_k(idx));
         %[p, pRank] = ssqls(sparse(Z'*B_k(idx,idx)*Z), -Z'*g_k(idx));
         [p, pRank]  = sls(sparse(Z'*B_k(idx,idx)*Z), -Z'*g_k(idx));
      end
   else
      if isempty(Z)
         [Q, R, E, pRank] = ComputeQR(B_k(idx,idx), epsRank);
      else
         [Q, R, E, pRank] = ComputeQR(Z'*B_k(idx,idx)*Z, epsRank);
      end

      if pRank==0
         pCond=Inf;
      else
         pCond=abs(R(1,1))/max(1E-200,abs(R(pRank,pRank)));
      end
      if pRank == 0
         if any(isnan(g_k)) | any(any(isnan(B_k)))
            fprintf('ERROR!!! THE USER HAS RETURNED NAN ELEMENTS ')
            fprintf('IN THE RESIDUAL OR JACOBIAN\n')
         end
      end

      if isempty(Z)
         p = tomsol(6, Q, R, E, pRank, -g_k(idx));
      else
         p = full(Z * tomsol(6, Q, R, E, pRank, full(-Z'*g_k(idx))));
      end
   end

elseif Method == 1 & GN % Gauss-Newton step with Singular Value Decomposition
   if isempty(Z)
      [U S V] = svd(full(J_k(:,idx)));
      k = length(idx);
   else
      [U S V] = svd(full(J_k(:,idx)*Z));
      k = size(Z,2);
   end
   S_inv = zeros(size(S,2),size(S,1));
   S_inv(1,1) = 1/S(1,1);
   pRank=1;
   for i = 2:k
       if S(i,i) > epsRank*S(1,1)
          S_inv(i,i) = 1/S(i,i);
          pRank=i;
       end
   end
   if S(pRank,pRank)==0
      pCond=Inf;
   else
      pCond=S(1,1)/max(1E-200,S(pRank,pRank));
   end
   if isempty(Z)
      p = V * (S_inv * (U' * (-r_k)));
   else
      p = full(Z * V * (S_inv * (U' * (-r_k))));
   end

%fprintf('Projected gradient %25.14f\n',Z'*g_k(idx));
%fprintf('Directed derivative %25.14f\n',g_k(idx)'*p);
%xprint(p,'p');

elseif Method == 1 & ~GN % In Fletcher-Xu hybrid method. Gauss-Newton-BFGS
   %-------------------------------------------
   % Solve normal equations   B_k*p=-g_k. if GN=0 then B_k is BFGS-update
   if isempty(Z)
      [U S V] = svd(full(B_k(idx,idx)));
      k = length(idx);
   else
      [U S V] = svd(full(Z'*B_k(idx,idx)*Z));
      k = size(Z,2);
   end
   S_inv = sparse(size(S,2),size(S,1));
   S_inv(1,1) = 1/S(1,1);
   pRank=1;
   for i = 2:k
       if S(i,i) > epsRank*S(1,1)
          S_inv(i,i) = 1/S(i,i);
          pRank=i;
       end
   end
   if S(pRank,pRank)==0
      pCond=Inf;
   else
      pCond=S(1,1)/max(1E-200,S(pRank,pRank));
   end
   if isempty(Z)
      p = full(V * (S_inv * (U' * (-g_k(idx)))));
   else
      p = full(Z * V * (S_inv * (U' * (-Z'*g_k(idx)))));
   end
elseif Method == 2 & GN
   % Gauss-Newton with Matlabs inversion routine (Uses QR)
   if isempty(Z)
      p =  -J_k(:,idx) \ r_k;
      pRank = rank(full(J_k(:,idx)));
   else
      p = full(Z*( (J_k(:,idx)*Z) \ (-r_k)));
      pRank = rank(full(J_k(:,idx)*Z));
   end
elseif Method == 2 & ~GN
   % Hybrid method. Solving normal equations;
   % Using Matlabs inversion routine (Uses QR)
   if isempty(Z)
      p =  -B_k(idx,idx) \ g_k(idx);
      pRank = rank(full(B_k(idx,idx)));
   else
      p = full( (Z'*B_k(idx,idx)*Z) \ (-Z'*g_k(idx)));
      pRank = rank(full(Z'*B_k(idx,idx)*Z));
   end
elseif Method == 3 & GN
   % Gauss-N; Explicit computation of pseudoinverse, pinv(J_k)
   if isempty(Z)
      p = pinv(full(J_k(:,idx)),epsRank)*(-r_k);
      pRank = rank(full(J_k(:,idx)));
   else
      p = full(Z*(pinv(full(J_k(:,idx)*Z),epsRank)*(-r_k)));
      pRank = rank(full(J_k(:,idx)*Z));
   end
elseif Method == 3 & ~GN
   % Hybrid method. Solving normal equations;
   % Explicit computation of pseudoinverse, pinv(J_k)
   if isempty(Z)
      p = pinv(full(B_k(idx,idx)),epsRank)*(-g_k(idx));
      pRank = rank(full(B_k(idx,idx)));
   else
      p = full((pinv(full(Z'*B_k(idx,idx)*Z),epsRank))*(-Z'*g_k(idx)));
      pRank = rank(full(Z'*B_k(idx,idx)*Z));
   end
end

% =========================================
function Prob = CheckSep(Prob)
% =========================================
% Check if separable nonlinear least squares
if isfield(Prob.LS,'SepAlg')
   if Prob.LS.SepAlg
      parr=nargout(Prob.FUNCS.r);
      if parr < 4
         fprintf('clsSolve: ERROR!!! Flag Prob.LS.SepAlg=1, but ');
         fprintf('%s',Prob.FUNCS.r);
         fprintf(' has only %d output parameters\n',parr);
         fprintf('Must return r, J, z and Jz, at least 4 parameters!\n');
         fprintf('Try to run ordinary nonlinear least squares.\n');
         Prob.LS.SepAlg=0;
      end
   end
end


% =========================================
function ConTest(b, A, x_k, me, bTol, PriLev)
% =========================================

% check if x_k is feasible (might be skipped)
test = (b-A*x_k) > bTol;
contest(1:me)=0;
if any(contest)
   if PriLev >=0
      fprintf('\n ERROR ! x_k not feasible')
      fprintf(', Violated constraint:\n')
      xprinti(find(contest),'index:');
      xprint(x_k,' x_k ');
   end
end

% =========================================
function [x_k, me, mc, mi, LCflag, A, b, W_kset]=leqCheck(Prob, ...
         SolverQP, x_k, b_L, b_U, x_L, x_U, bEqual, Ax, bTol, PriLev)
% =========================================

LCflag = 0;

eqCon=find(bEqual & ~isinf(b_L));
ineqLow=find(~bEqual & ~isinf(b_L));
ineqUpp=find(~bEqual & ~isinf(b_U));

me = sum(bEqual); % Number of equality constraints
mc = length(eqCon)+length(ineqLow)+length(ineqUpp);
mi = mc-me; % Number of inequality constraints

if mc > 0 % If any linear constraints

   LCflag = 1;

   % Check if x_0 is feasible due to the linear constraints.
   % If not, minimize the sum of squares of relative deviation
   % between x(i) and x_0(i) by solving:
   %
   %    min  0.5*(x-x_0)'*B*(x-x_0)           min   0.5*x'*B*x - x_0'*B*x
   %    s/t   x_L <=   x  <= x_U       <=>    s/t   x_L <=   x  <= x_U
   %          b_L <=  Ax  <= b_U                    b_L <=  Ax  <= b_U

   if any( Ax + bTol < b_L | Ax - bTol > b_U ) % Not feasible
      n = length(x_k);
      % Should we use relative weighting or not.
      % May cause ill-conditioning if to high F-values
      if n > 1000
         F = spdiags(1./max(1E-3,x_k.^2),0,n,n);
      else
         F = diag(1./max(1E-3,x_k.^2));
      end
      % Setup structure used in QP call
      ProbQP = qpAssign(F,-F'*x_k,Prob.A,b_L,b_U,x_L,x_U,zeros(n,1),'QP-Feas');
      ProbQP.optParam = optParamDef(SolverQP,2,n,0, size(Prob.A,1));
      ProbQP.optParam.MaxIter = max([500,3*n,ProbQP.optParam.MaxIter]);

      ResultQP = tomRunMini(SolverQP,ProbQP);

      x_k = ResultQP.x_k;

      if PriLev > 0
         fprintf('\n Given x_0 not feasible due to Linear Constraints.\n');
         xprint(x_k,' New x_0:',' %14.14e');
      end
   end

   % Put linear equalities first in A and turn '<=' inequalities
   % to '>=' inequalities

   A = [Prob.A(eqCon,:);Prob.A(ineqLow,:);-Prob.A(ineqUpp,:)];
   b = [b_L(eqCon);b_L(ineqLow);-b_U(ineqUpp)];

   W_kset = [ones(me,1);zeros(mi,1)];
else
   A=[]; b=[];
   W_kset = [];
end

% ==========================================================
function [set, nract]=x_kActive(nract0, x_k, x_L, x_U, xTol, Arank)
% ==========================================================


n=length(x_k);
% This is the old variant, not well thought.
%Free=n-Arank-nract0;

Free=n-Arank;
set = zeros(n,1);
set(abs(x_k-x_L) <= xTol * max(1,abs(x_k)))=-1;

% Will now set to 1 when x_L == x_U
set(abs(x_k-x_U) <= xTol * max(1,abs(x_k)))=1;
nract = sum(set~=0);

if nract > Free
   % Must reduce the number of active variables.
   % This simple strategy might lead to reduction from full rank
   ix = find(set~=0 & x_L~=x_U);
   set(ix(1:nract-Free))=0;
   nract=sum(set~=0);
end

%set = x_L==x_U;
%for i=1:n
%    if ~set(i)
%       if x_k(i) >= x_U(i)	% x active on upper bound
%          x_k(i) = x_U(i);      % Fix x exactly at bound
%          if nract < n-Arank
%             set(i) = 1;
%             nract = nract + 1;
%          end
%       elseif x_k(i) <= x_L(i)	% x active on lower bound
%          x_k(i) = x_L(i);      % Fix x exactly at bound
%          if nract < n-Arank
%             set(i) = -1;
%             nract = nract + 1;
%          end
%       else
%          if x_k(i) < x_L(i)+xTol	  % x close to lower bound
%             x_k(i) = x_L(i);             % Move x to bound
%             if nract < n-Arank
%                set(i) = -1;
%                nract = nract + 1;
%             end
%             changed = 1;
%          elseif x_k(i) > x_U(i)-xTol     % x close to upper bound
%             x_k(i) = x_U(i);             % Move x to bound
%             if nract < n-Arank
%                set(i) = 1;
%                nract = nract + 1;
%             end
%             changed = 1;
%          else
%             set(i) = 0;
%          end
%       end
%    end
%end
%nract = sum (set~=0); % Number of active variables, i.e. on bound.

% ==========================================================
function itPrint(PriLev, x_k, f_k, g_k, r_k, J_k, k, pRank, pCond, nract, ...
                 set, zSum, KT, KTT)
% ==========================================================

if PriLev > 1
   fprintf('==================================================\n');
   fprintf('Iteration no: %4.0f  Func %30.20f pRank %d',k,f_k,pRank);
   if ~isnan(pCond), fprintf(' pCond %10.5e',pCond); end
   n=length(x_k);
   if n < 50,
      fprintf(' Cond %10.5e',cond(full(J_k)));
   elseif n > 50
      if issparse(J_k)
         %fprintf(' Cond %10.5e',condest(J_k));
      else
         %fprintf(' Cond %10.5e',condest(sparse(J_k)));
      end
   end
   fprintf('\n');
   fprintf('Sum of infeasibilities             %30.20f\n',zSum);
   if ~(isinf(KT) | isnan(KT))
      fprintf('Kuhn-Tucker point %d. K-T Error     %30.20f\n',KTT,KT);
   end
end
if PriLev > 2
   xprint(x_k,'x_k:');
   if nract > 0
      fprintf('Active variables = %3.0f\n',nract);
      xprint(set,'Set:',' %2.0f',25);
   end
   xprinte(g_k,'g_k:');
   if PriLev > 3, xprint(r_k,'r_k:'); end
   if PriLev > 4, mPrint(J_k,'J_k',' %13.6e',5); end
end

% ==========================================================
function A=secUpdat(s,y,B,v)
% ==========================================================
%
% secUpdat computes the least change formulation of secant updates.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 27, 1997. Last modified June 22, 1999.
%

vTs=v'*s;
z=(y-B*s)/vTs;
A=z*v' + v*z' - (z'*s/vTs)*v*v';

% ==========================================================
function [h_k,cErr] = ConstraintError(v1,v2,alg)
% ==========================================================
% Only called from function KuhnTucker

if isempty(v1) & isempty(v2)
   cErr=zeros(0,1);
   h_k=0;
   return
end

cErr=max(v1,v2);

if alg < 2
   h_k = norm(max(0,cErr),Inf);  % Constr. violation, max norm
else
   h_k = norm(max(0,cErr),1);    % Constr. violation, L1 - sum of abs
end
% ==========================================================
function [KT, KTT, h_k, cErr] = KuhnTucker(Prob, Result,...
    bTol, zTol, xTol, x_k, x_L, x_U, g_k, z, bl, bu, xEqual, zEqual)
% ==========================================================

[v_k, Zv, P, Z_L, cErr, ceq, cineq, gProj0] = LagMult(Prob,Result);

n=length(x_k);
mTot=length(z);
if mTot == 1
   Anorm=norm(Prob.A);
elseif mTot==0
   Anorm=[];
else
   Anorm=sqrt(sum(Prob.A'.^2));
end
if mTot > 0
   myMax=max([norm(g_k),max(abs(v_k(1:n))) max(Anorm(:).*v_k(n+1:n+mTot))]);
   if myMax==0, myMax=1; end

   KT=norm(g_k-v_k(1:n)-Prob.A'*v_k(n+1:n+mTot))/myMax;
else
   myMax=norm(g_k);
   if myMax==0, myMax=1; end
   KT=(g_k-v_k(1:n))/myMax;
end
KT = max(KT);

%fprintf('Kuhn-Tucker error %30.20f\n',KT);

[h_k cErr] = ConstraintError(z-bu,bl-z, 0);

if KT <= bTol & (mTot==0 | all(cErr <= zTol.*max(1,abs(z))))
   xT=xTol*max(1,abs(x_k));
   iL=abs(x_L-x_k) < xT & ~xEqual;
   iU=abs(x_U-x_k) < xT & ~xEqual;
   iF=find(~(iL | iU | xEqual));
   KTT=all(v_k(iL) >= -xT(iL)) & all(v_k(iU) <= xT(iU)) & ...
       all(abs(v_k(iF)) < xT(iF));
%if KTT==0
%KTT
%iL
%iU
%iF
%   all(v_k(iL) >= -xT(iL))
%all(v_k(iU) <= xT(iU))
%       all(abs(v_k(iF)) < xT(iF))
%[x_L x_k x_U v_k(1:n)]
%pause
%end
   if KTT & mTot > 0
      zT=zTol.*max(1,abs(z));
      iL=abs(bl-z) < zT & ~zEqual;
      iU=abs(bu-z) < zT & ~zEqual;
      iF=find(~(iL | iU | zEqual));
      KTT=all(v_k(n+find(iL)) > -zT(find(iL))) & ...
          all(v_k(n+find(iU)) <  zT(find(iU))) & ...
          all(abs(v_k(n+iF)) < zT(iF));
   end
else
   KTT=0;
end
% ==========================================================
function LineResult=Armijo(x_0,p,f_0,g_0,rho,sigma,alphaMax,Prob,varargin)
% ==========================================================
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 27, 1997. Last modified Feb 15, 2005.
%

%MaxIter = min(4,1+ceil(log(1E-14)/log(rho)));
MaxIter = 1+ceil(log(1E-14)/log(rho));
fp_0    = g_0'*p;
    %r_k = nlp_r(x_0, Prob, varargin{:} );    % residual value
    %J_k = nlp_J(x_0, Prob, varargin{:} );    % residual value
    %g_k=J_k'*r_k
alpha     = min(1,alphaMax);
alphaV    = ones(MaxIter+1,1);
alphaV(1) = alpha;
    % NOT USED f_0S    = sqrt(2*f_0);
alg = Prob.Solver.Alg;
if alg==7
   LineResult.J_k = [];
   LineResult.g_alpha = g_0;
end

%fp_0
for i = 1:MaxIter
    x_k = x_0 + alpha * p;
    r_k = nlp_r(x_k, Prob, varargin{:} );     % residual value
    f_k = 0.5*(r_k'*r_k);
    %f_k = nlp_f(x_k, Prob, varargin{:} );     % function value
    %format long
    %disp([i, f_k, f_0 + sigma*alpha*fp_0,f_k-f_0 + sigma*alpha*fp_0])
    %if f_k <= f_0 + sigma * alpha * fp_0
    if f_k-(f_0+sigma * alpha * fp_0) <= max(1,abs(f_0))*1E-6
       LineResult.alpha = alpha;
       LineResult.f_alpha = f_k;
       %r_k = nlp_r(x_k, Prob, varargin{:} );     % residual value
       LineResult.r_k = r_k;
       if alg ~= 7
          J_k = nlp_J(x_k, Prob, varargin{:} );  % Jacobian value
          LineResult.J_k = J_k;
          if ~isempty(J_k)
             LineResult.g_alpha = J_k'*r_k;
          else
             LineResult.g_alpha = [];
          end
       end
       LineResult.alphaVec = alphaV(1:i);
       return
    end
    alpha = alpha*rho;
    alphaV(i+1) = alpha;
end

LineResult.alpha = alpha;
LineResult.f_alpha = f_k;
%r_k = nlp_r(x_k, Prob, varargin{:} );     % residual value
LineResult.r_k = r_k;
if alg ~= 7
   J_k = nlp_J(x_k, Prob, varargin{:} );  % Jacobian value
   LineResult.J_k = J_k;
   if ~isempty(J_k)
      LineResult.g_alpha = J_k'*r_k;
   else
      LineResult.g_alpha = [];
   end
end
LineResult.alphaVec = alphaV(1:MaxIter);

% ==========================================================
function [x_k, pRank]  = sls(C, y, damp, aTol, bTol, condLim, MaxIter, x_0)
% ==========================================================

if nargin < 8
   x_0 = [];
   if nargin < 7
      MaxIter = [];
      if nargin < 6
         condLim = [];
         if nargin < 5
            bTol = [];
            if nargin < 4
               aTol = [];
               if nargin < 3
                  damp = [];
end, end, end, end, end, end

PriLev=0;

% Linear least squares

m = length(y);
n = size(C,2);

if isempty(damp), damp = 0; end

%aTol       = Prob.optParam.xTol;
%bTol       = Prob.optParam.bTol;
%MaxIter    = 100*(n+k);
%MaxIter    = Prob.optParam.MaxIter;

WantStdErr = 1;
D          = [];

% Can set nOut to a filename, or Fortran file unit (writes to lsqrout.txt)
nOut = [];

[ x_k, Inform, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
    Tlsqr( m, n, C, [], [], y, damp, aTol, bTol, condLim, ...
    MaxIter, WantStdErr, nOut, D, x_0 );

pRank = sum(abs(x_k) > 1E-14);

% Linear least squares
%r   = C*x_k-y;
%f_k = 0.5*(r'*r);
%g_k = C'*r;

if PriLev > 0

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

ExitText = Text;

   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLSQR solving QR Problem \n');
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LSQR: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   fprintf('Sum of squares at solution x %26.18f\n\n',f_k);
   fprintf('Iterations %11d\n',Iter);
   fprintf('Variables  %11d\n',n);
   % fprintf('Slacks     %11d\n',k);
   % fprintf('Equations  %11d\n',size(A,1));
   fprintf('rNorm   %13.7e\n',rNorm);
   fprintf('xNorm   %13.7e\n',xNorm);
   fprintf('aNorm   %13.7e\n',aNorm);
   fprintf('aCond   %13.7e\n',aCond);
   fprintf('arNorm  %13.7e\n',arNorm);
   fprintf('\n');
end

% ==========================================================
function [x] = ItRef(A, b, x, R, pRank, MaxIter)
% ==========================================================

for i=1:MaxIter
    r   = b-A(:,1:pRank)*x(1:pRank);
    rNorm = norm(r);
    %ATb = A'*r;
    %y   = R'\ATb;
    y   = tomsol(2,R',A(:,1:pRank)'*r);
    %dx  = R\y;
    dx  = tomsol(3,R,y);
    x(1:pRank)=x(1:pRank)+dx;
end

r   = b-A(:,1:pRank)*x(1:pRank);
NewNorm = norm(r);

% MODIFICATION LOG:
%
% 980409  mbk   Gradient convergence test on projected gradient.
%               Reactivation tolerance for alpha changed to eps_con.
%               Null space update made before convergence tests instead of after.
%               Residual and Jacobian computation on line ~ 261 aborted.
% 980413  mbk   Temporary debug output
%
% 980429  hkh   Changed to clsSolve and made modifications for new TOMLAB,
%               using structure prob.
%
% 980505  mbk   Using logic vector fixbnd together with lambda to determine if
%               performing check for convergence.
%               Convergence check "norm(x_k-x_km1) <= eps_x" changed to
%               "( norm(x_k-x_km1) <= eps_x )&(alphaMax > 1E-14)" to avoid
%               flagging for convergence when iteration points are close because
%               of restriction on alphaMax.
%
% 980506  hkh   Avoid global r_k and J_k. Instead extra pars from line search
%               Added initialization W_kset = []; for unconstrained problems
%
% 980514  hkh   result renamed to Result. Result initialized by ResultDef.
%               prob   renamed to Prob.
%               Using xState, bState like NPSOL etc.
% 980515  hkh   Define which solver in structure: Result.Solver='clsSolve';
%
% 980608  mbk   Modified for linear constraints
% 980608  mbk   Computation of mc corrected.
%
% 980608  hkh   r_k, J_k not always computed in linesrch. Check and update.
%
% 980608  mbk   Call to presolve analysis routine.
% 980816  hkh   Changed qpiold to qpiOld
% 980922  hkh   Change from optPar to structure optParam. Call LineSearch.
% 981005  mbk   Changes in comments concerning type of convergence.
% 981005  mbk   clsDef instead of ucDef.
% 981005  hkh   New call to LineSearch. Using LinePlot. Cleanup.
%                New flags linePlot and DEBUG.
% 981006  mbk   alphaMax set to 1E20 in init for conv. check in 1:st it.
% 981013  hkh   Added call to iniSolve and endSolve
% 981017  hkh   Deleted the Jacobian as argument to LinePlot
% 981026  hkh   Solver is name of solver. SolverAlgorithm is description
% 981027  hkh   Change printing levels. Use Prob.Solver.Alg to define alg
% 981028  hkh   Put f_0 = f(x_0) as field in Result. Test on SepAlg true.
% 981102  hkh   Check if field SepAlg is defined before checking on it.
% 981105  hkh   Unnecessary to set Result.Prob here, done in endSolve
% 981108  hkh   Check convergence always if fRed0cnt > LowIts, otherwise
%               infinite loop (to MaxIter). Proj.gradient may be low.
% 981108  hkh   Add fourth parameter z in call when SepAlg true
% 981110  mbk   Don't care about negative lambdas corresponding to fixed
%               variables when determine which constraint or which variable
%               to release. Use lambda_hat.
%               Change test of fRed0cnt to use 100 * machine precision
% 981112  mbk   eps_con replaced by cTol. eps_con was set to cTol before.
%               Call to preSolve moved to the beginning of the code.
%               qpSolve is called before converting to Ax=b,Ax>=b format.
% 981120  mbk   b_k=JJ; changed to B_k=JJ; if alg==0.
%               B_k is now initially set to eye(n) or optParam.QN_InitMatrix
%               not to J'J as before.
% 981128  hkh   Change printing levels for some output
% 981130  hkh   Use hard coded flag if to use QPOPT or qpSolve
% 981208  hkh   Conflict between A=Prob.A and expanded A in Ax >= b.
% 981208  hkh   All Lagrange multipliers not returned. Now setting
%               Result.v_k=[lam_L(:);lambda(:)];
% 990216  hkh   Check for rank(A)==0, due to variables on bounds.
% 990221  hkh   Print p_rank, compute pCond and print.
% 990221  hkh   Check if Prob.FUNCS.d2r is defined, before computing 2nd der
% 990330  hkh   Compute exact pseudo rank pRank as output, and safe pCond
% 990910  hkh   Add IterPrint one row printing each iteration
% 000909  hkh   Use optParam.bTol, not cTol, as tol for linear constraints
% 000910  hkh   Skip using optParam.NOT_release_all, hard coded instead as 0
% 000916  hkh   Adding text for ExitText
% 000923  hkh   Clean up. Divide between optParam and LineSearch
%               Move fLow to fLowBnd in LineSearch, LineAlg to LineSearch
% 000927  hkh   Change to Prob.PriLevOpt
% 001107  hkh   Use only bTol
% 011205  hkh   Add Inform entry 11
% 020409  hkh   Use Prob.nState and Prob.Mode
% 020417  hkh   Add the use of sqr2 to efficiently store the Q matrix
%               New routine ssqls is called.
% 020531  hkh   LargeScale = 0, if n<10, to avoid bugs in sqr2 package
% 030510  hkh   Add LineParam.LineAlg parameter, LineParam.agFac
% 030510  hkh   Add LineParam.agFac, LineParam.sigma for Armijo-Goldstein
% 030510  hkh   Add output Result.IterGN, sum of Gauss-Newton steps
% 030515  hkh   Bug in Huschens method, A_k was returned after update
% 030515  hkh   Addition of Structured MBFGS method by Wang-Li-Qi
% 030924  hkh   Use spdiags and not diag if > 1000 variables to formulate QP
% 040111  hkh   Change call to inisolve
% 040414  hkh   Bug when calling ssqls with J_k(:,idx)*Z, Z*p then needed
% 040602  med   PriLevOpt added
% 040604  hkh   Make Z and B_k sparse
% 040917  hkh   Avoid gProj test if gProj and idx empty. ; missing if Inform=32
% 040920  frhe  Logicals casted to doubles to prevent mult's with logicals.
% 041203  hkh   Serious bugs in Method 2 and 3 computing search direction
% 041207  hkh   Print: Convergence 5 if PriLev >1, not >-1.
% 041208  hkh   Change two PriLev >= 0 to > 0
% 050117  med   mlint review
% 050121  frhe  Argument order in nlp_d2r calls changed.
% 050205  hkh   ExitFlag = 0, not 1, for eps_absf convergence
% 050221  hkh   Check if NaN or Inf in r_k at x_0, then exit safely
% 050221  hkh   Check if NaN or Inf in J_k at x_0, then exit safely
% 050221  hkh   Revised Wang, Li, Qi Structured MBFGS method
% 050221  hkh   Made Li-Fukushima MBFGS method, slight change from Wang,Li,Qi
% 050304  hkh   if f_red < 1E-30*max(1,f_k), try sigma 0.1 & 0.01, before stop
% 050727  hkh   Added Broydons method
% 051216  med   Extra inform text for 18 removed, help updated
% 060216  hkh   If pRank==0, check for NaN and display Error message
% 060804  hkh   Skip ssqls, sqr2 sparse QR package, non-working for new Matlab
% 060804  hkh   For LargeScale=1, use subfunction sls, calling Tlsqr
% 060804  hkh   Add iterative refinement routine ItRef, not used right now
% 060814  med   FUNCS used for callbacks instead
% 060814  hkh   Diff r_k-r_km1, not f_k-f_km1 in Broydens method,alg=7
% 060814  hkh   r_k not updated correctly in Broydens method
% 060814  hkh   Armijo linesearch revised, use alphaMax, no return of J_k,alg=7
% 060817  hkh   Always sparse QR when LargeScale, if LargeScale > 0, Method = 0;
% 060817  hkh   Return after error messages if any nonlinear constraints
% 060818  hkh   Use Prob.f_Low, changed comments about f_Low, fLowBnd
% 060818  hkh   Use Prob.f_Low instead of optParam.eps_absf for target f test
% 060818  hkh   Use optParam.eps_absf for close to target test
% 060818  hkh   Modified Convergence 3 test and Stop 2
% 070221  med   Help updated
% 070725  frhe  Estimate Jacobian when using Broyden's method.
% 070907  hkh   SolverQP picked from list, 1st with license, avoid GetSolver
% 080607  hkh   Use tomRun, not tomSolve. Move QP to leqCheck, no CreateProbQP
% 080607  med   Switched to tomRunMini
% 080619  hkh   Set p=zeros(length(idx),1) if pRank == 0, not n, avoid crash
% 090717  med   Residual calculation updated
% 110722  hkh   Added comments about weighting the NLLS problem
