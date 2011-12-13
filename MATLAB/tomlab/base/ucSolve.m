% ucSolve.m
%
% UCSOLVE is a solver implementing several algorithms for
% unconstrained minimization with bound constraints,
%
%	min      f(x)
%
%            x_L <= x <= x_U
%
% function Result = ucSolve(Prob, varargin)
%
% UCSOLVE implements the following algorithms:
%
% Algorithm: 0 = Default algorithm (BFGS or Newton)
%    1 = Newton with subspace minimization, using SVD
%    2 = Safeguarded BFGS with standard inverse Hessian update
%    3 = Safeguarded BFGS with Hessian update, and SVD or LU to solve
%    4 = Safeguarded DFP  with standard inverse Hessian update
%    5 = Safeguarded DFP  with Hessian update, and SVD or LU to solve
%    6 = Fletcher-Reeves CG
%    7 = Polak-Ribiere CG
%    8 = Fletcher conjugate descent CG-method
%
% Set the field Prob.Solver.Alg to one of the numbers above
%
% Note: The accuracy in the line search is critical for the performance
% of quasi-Newton BFGS and DFP methods, as well as the CG methods.
% If the accuracy parameter, Prob.LineParam.sigma == 0.9 (standard default),
% it is changed by ucSolve to:
%       Prob.Solver.Alg == 4,5   => sigma =  0.2
%       Prob.Solver.Alg == 6,7,8 => sigma =  0.01
%
% Bound constraints treated as described in
% Gill, Murray, Wright: Practical Optimization, Academic Press, 1981.
%
% Solver.Method    Method to solve equation system (Solver.Alg in [0,5])
%   0:  SVD (default).
%   1:  LU-decomposition.
%   2:  LU-decomposition with pivoting.
%   3:  Matlab built in QR. (not recommended)
%   4:  Matlab inversion. (not recommended)
%   5:  Explicit inverse. (not recommended)
% Solver.Method    Restart or not for C-G method: (Solver.Alg in [6,8])
%   0:  No restart in CG-method each n:th step (default).
%   1:  Use restart in CG-method each n:th step.
%
% INPUT PARAMETERS
%
%   Use conAssign.m (or probAssign) to initialize the Prob structure in the
%   TOMLAB format.  Fields used in structure Prob:
%
%    FUNCS.f   The routine to compute the function value, given as a string
%    FUNCS.g   The routine to compute the gradient vector, given as a string
%              If empty, or Prob.NumDiff ~=0 numerical differences are used
%              If Prob.ADObj = 1,  MAD automatic differentiation is used
%    FUNCS.H   The routine to compute the Hessian, given as a string
%              Only used for Newtons method.
%              If empty, or Prob.NumDiff < 0 numerical differences are used
%              If Prob.ADObj = -1, MAD automatic differentiation is used
%    x_0       Starting point
%    x_L       Lower bounds for x
%    x_U       Upper bounds for x
% PriLevOpt    Print level
%
%    Solver.Alg    Algorithm number
%    Solver.Method Search direction solution technique
%
%    optParam   Structure with optimization parameters
%       eps_f     Relative change in f, convergence tolerance
%                 ucSolve will stop if the test is true for LowIts iterations
%                 LowIts = 10; now, to change it, change on line 182.
%       eps_g     Gradient convergence tolerance
%       eps_Rank  Rank test tolerance
%       eps_x     Convergence tolerance in x
%       MaxIter   Maximal number of iterations
%       optParam.QN_InitMatrix  Initial Quasi-Newton matrix, empty ==> eye(n),
%       IterPrint Print short information each iteration
%       size_x    Approximate size of optimal variable values, normally 1
%       size_f    Approximate size of optimal function value, normally 1
%       xTol      Tolerance to judge if x is close to bounds
%                               otherwise use identity matrix
%     LineParam   Line search parameters, see LineSearch.m
%
%  PriLev:   Print level:  0 None, 1 Final result, 2 Each iteration, short
%            3 Each iteration, more info, 4 Line search info and Hessian
%
% The Prob structure could be created in the TOMLAB Quick format with
% calls to conAssign.m (or probAssign and tomFiles),
% or in the Init File format, see Users Guide.
%
% ---------------------------------------- extra parameters:
% VARARGIN: User defined parameters passed to f,g,H, r and J, and c and dc.
%
% OUTPUT PARAMETERS
%
% Result      Structure with results from optimization
%    x_k      Optimal point
%    v_k      Lagrange multipliers for bound constraints
%    f_k      Function value at optimum
%    g_k      Gradient vector at optimum
%    x_0      Starting value vector
%    H_k      The Hessian matrix, if computed (Newtons method)
%    B_k      The Quasi Newton matrix
%    xState   Variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%    Iter     Number of iterations
%    ExitFlag Flag giving exit status
%    ExitTest Text string giving ExitFlag and Inform information
%    Inform   Code telling type of convergence
%             1   Iteration points are close.
%             2   Projected gradient small.
%             3   Iteration points are close and projected gradient small.
%             4   Relative function value reduction low for 10 iterations.
%             5   Iteration points are close and relative function value
%                 reduction low for 10 iterations.
%             6   Projected gradient small and relative function value
%                 reduction low for 10 iterations.
%             7   Iteration points are close, projected gradient small and
%                 relative function value reduction low for 10 iterations.
%           101   Max no of iterations reached.
%           102   Function value below given estimate.
%           103   Convergence to saddle point (eigenvalues computed)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1992-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Dec 2, 1992.    Last modified Jul 20, 2009.

function Result = ucSolve(Prob, varargin)

if nargin < 1
   error('ucSolve needs input structure Prob');
end

global n_f xLast

solvType=checkType('uc');

Prob=ProbCheck(Prob,'ucSolve',solvType);

% Default is Newton or BFGS

Alg=Prob.Solver.Alg;
if isempty(Alg), Alg=0; end
if Alg == 0
   if isempty(Prob.NumDiff) | Prob.NumDiff==0
      if ~isempty(Prob.FUNCS.g)
         % Try with Newton, when analytic derivatives
         Alg=1;
      else
         Alg=2;
      end
   else
      Alg=2;
   end
end

if Alg == 1
   Prob = iniSolve(Prob,solvType,2,0);
else
   Prob = iniSolve(Prob,solvType,1,0);
end
% SICK - to avoid output in case misuse of ucSolve with nonlinear constraints
Prob.ConsDiff = 0;

DEBUG=0;     % Debugging flag
% If true (==1) plot line search problem at each iteration
plotLine=Prob.plotLine;

% Pick up input variables from Prob structure

% Gateway routine names (text strings)
f    = 'nlp_f';
g    = 'nlp_g';

H_k  = [];
B_k  = [];

[n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

lam_L = zeros(n,1);     % Initial zeros for Lagrange multipliers

optParam  = Prob.optParam;
LineParam = Prob.LineParam;
fLow      = Prob.f_Low;
% Lower bound on function value, not critical that it is tight
LineParam.fLowBnd = max(Prob.LineParam.fLowBnd,fLow); 

eps_x    = optParam.eps_x;    % Convergence tolerance in x
eps_g    = optParam.eps_g;    % Gradient convergence tolerance
epsRank  = optParam.eps_Rank; % Rank test tolerance
MaxIter  = optParam.MaxIter;  % Maximal number of iterations
size_x   = optParam.size_x;   % Approximate size of optimal variable values 
size_f   = optParam.size_f;   % Approximate size of optimal function value 

% If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bounds
xTol =  optParam.xTol;    
% If (f_km1-f_k) < eps_f * max(f_k,size_f) for LowIts == 10 iter, stop
eps_f = optParam.eps_f;     
% No of iterations with low reduction before convergence
LowIts = 10;

if LineParam.LineAlg > 1
   CURV=1;
   LineParam.LineAlg=LineParam.LineAlg-2;
   xLast=x_km1;
else
   CURV=0;   % No curvilinear search
end

Result=ResultDef(Prob);

method    = Prob.Solver.Method;   % Search direction solution technique

if isempty(method), method=0; end

PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration

sigma  = DefPar(LineParam,'sigma',0.9);
if sigma == 0.9
   if (Alg==4 || Alg==5) 
      % Must make more accurate line search for DFP methods
      LineParam.sigma = 0.2;
      if PriLev > -1
         fprintf('Reset line search parameter sigma from 0.9 to %7.3f\n', ...
                  LineParam.sigma);
      end
   end
   if any(Alg==[6 7 8]) 
      % Must make more accurate line search for C-G methods
      LineParam.sigma = 0.01;
      if PriLev > -1
         fprintf('Reset line search parameter sigma from 0.9 to %7.3f\n', ...
                  LineParam.sigma);
      end
   end
end

[SolvAlg, Alg, method] = NameAlg(Alg, method);

Result.Solver='ucSolve';
Result.SolverAlgorithm=SolvAlg;

if Alg==2 || Alg==4
   method=6; % Explicit Inverse Hessian. No equations to solve.
elseif any(Alg==[6 7 8 ])
   beta=0;
   cg_step=1;
   if method == 8 || method == 1
      cg_restart=1; % Use restart in CG-method each n:th step.
   else
      cg_restart=0; % No restart, DEFAULT
   end
   method=7;
else
   if method > 5, method=0; end % Safety check if wrong input
end

% Find initial Quasi-Newton matrix
if any(Alg==(2:5)) 
   if ~isempty(optParam.QN_InitMatrix)
      B_k = optParam.QN_InitMatrix;   % B = User given matrix
      if size(B_k,1)~=n || size(B_k,2)~=n 
         disp('Illegal dimensions on initial Quasi Newton matrix');
         disp('ucSolve is using identity matrix instead');

         B_k = eye(n);                % B = Use identity matrix
      end
   else
      B_k = eye(n); % B = Start Quasi-Newton with identity matrix
   end
end
%
%	Init
%
k = 0;   % Iteration index. Initial value comp is step 0.

set = zeros(n,1);
EigStepTried=0;

x_km1 = inf*ones(n,1);
if CURV~=0
   xLast=x_km1;
end

alpha0cnt=0; stop = 0; f_km1=Inf; fred0cnt=0;
eig_step=0;

if IterPrint
   fprintf('Iteration Function         f(x)           |step|     ');
   fprintf('line search   |gradient|\n');
   fprintf('           Count                                       ');
   fprintf('   step                   \n');
end
p_full=0;
alpha=0;

while 1
% Check if variables near bound and set up active variable indicator set
   changed = 0;
   set_0=set;
   for i=1:n
       if x_k(i) >= x_U(i)	% x active on upper bound
          x_k(i) = x_U(i);      % Fix x exactly at bound
          set(i) = 1;
       elseif x_k(i) <= x_L(i)	% x active on lower bound
          x_k(i) = x_L(i);      % Fix x exactly at bound
          set(i) = -1;
       else
          if x_k(i) < x_L(i)+100*eps	  % x close to lower bound
             x_k(i) = x_L(i);             % Move x to bound
             set(i) = -1;
             changed = 1;
          elseif x_k(i) > x_U(i)-100*eps  % x close to upper bound
             x_k(i) = x_U(i);             % Move x to bound
             set(i) = 1;
             changed = 1;
          else
             set(i) = 0;
          end
       end
   end
   if k==0,set_0=set;end % 1st iteration; release from bounds directly
   nract = sum (set~=0); % Number of active variables, i.e. on bound.
   b_idx = find(set);    % Index vector for active variables, the fixed vars
   idx = find(~set);     % Index vector for nonactive variables, the free vars
   if changed || k==0
      Prob.nState = double(k==0);
      Prob.Mode   = 2;
      f_k = nlp_f(x_k,Prob, varargin{:});
      Prob.Mode   = 1;
      g_k = nlp_g(x_k,Prob, varargin{:});
      if isempty(g_k)
         Prob.NumDiff=1;
         g_k = nlp_g(x_k,Prob, varargin{:});
      end
      g_km1 = g_k;
      if k==0
         Result.x_0=x_k;
         Result.f_0=f_k;
      end
      Prob.nState = 0;
   end
   if PriLev > 1
      fprintf('==================================================\n');
      fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);
   end
   
   if IterPrint
      fprintf(' %5d  %7d   %20.17f %11.7e %10.6f %17.7e\n', ...
              k, n_f, f_k, norm(p_full), alpha, norm(g_k));
   end

   if PriLev > 2
      xprint(x_k,'x_k:');
      xprint(set,'set:',' %2.0f',25);
      xprinte(g_k,'g_k:');

   end

%
% Compute the Hessian, if Newtons method
%
   if Alg==1
      H_k = nlp_H(x_k,Prob, varargin{:});
      if isempty(H_k) && Alg==1
         % Must switch to BFGS from Newton
         disp('ucSolve: SWITCH TO BFGS, NEWTON IMPOSSIBLE')
         Alg=2;
         Result.SolverAlgorithm='Safeguarded BFGS';
         B_k = eye(n); % B = Start Quasi-Newton with identity matrix
      end
      if PriLev > 2 && (~isempty(H_k) || Prob.NumDiff==0)
         fprintf('Eigenvalues for Hessian H_k\n');
         xprinte(eig(full(H_k)),'eig:');
      end
      if PriLev > 3 & ~isempty(H_k)
         fprintf('The Hessian H_k\n');
         PrintMatrix(H_k,'H_k:','%10.6e',6);
      end
   end

%
% Check if any active variable shall be released
%
   xR=[]; xRp=[];  % Var# released and gradient value
   release=0;
   % Dangerous to use 1st order estimate. May ==> infinite loop
   if nract > 0 && alpha0cnt < 3
      lam_L = -set .* g_k; % Negative grad on upper bound, positive on lower
      if PriLev > 2
         fprintf('1st order Lagrange multiplier estimate. ');
         xprinti(b_idx,'Var#:');
         xprinte(lam_L(b_idx),'Lam:');
      end
      for i=1:nract
          j=b_idx(i);
          if nract < n
             % If lam_L > 0 optimal. Change sign in test
             % Release vars with negative Lagrange multiplier if they not
             % just hit a bound.
             % Also try a gradient step if releasing a vars result in 0 step.

             if (double((set_0(j)==set(j)) & (x_U(j)~=x_L(j)))*lam_L(j) < -xTol) |...
                ((set_0(j)==0) & double((alpha==0) & (x_U(j)~=x_L(j)))*lam_L(j) < -xTol)
                xR=[xR;j];
                xRp=[xRp;-g_k(j)];
                set(j)=0;
                release=1;
             %elseif (set_0(j)~=set(j))*(x_U(j)~=x_L(j))*lam_L(j) < -xTol
             %   % When parameter just reached a bound and signals
             %   % release, do not test for convergence
             %   %testconv=0;
             %   %disp('NOW CHANGED TO TESTCONV')
             %   if PriLev > 2
             %      disp('ucSolve: RELEASE DIRECTLY FROM BOUND')
             %   end
             %   xR=[xR;j];
             %   xRp=[xRp;-g_k(j)];
             %   set(j)=0;
             %   release=1;
             end
          else % No freedom to move. Must check directly if optimum
             if double((x_U(j)~=x_L(j)))*lam_L(j) < -xTol
                set(j)=0;
                release=1;
             end
          end
      end
   end

%
% Check convergence conditions
%

   Inform=0;

   if release 
      nract = sum (set~=0);
      b_idx = find(set);
      idx = find(~set); % Index vector for non active variables, the free vars
      if PriLev > 2
         fprintf('Release variable from bound. '); 
         xprinti(xR,'var#:');
      end

   else % Check convergence criteria

      %if all(lam_L >=-xTol) % Only check if possible optimum

      if max(abs(x_k-x_km1)./max(abs(x_k),size_x)) <= eps_x   
         if PriLev >= 1
            disp('*** Convergence 1, Iteration points are close ***');
         end
         Inform=Inform+1;
         flag=0;
         stop = 1;
      end
        
      if max(abs(g_k(idx)).* max(abs(x_k(idx)),size_x)) <= ...
           eps_g * max(abs(f_k),size_f)
         if PriLev >= 1
            % The projection is trivial in the case of simple bounds.
            disp('*** Convergence 2, Projected gradient small ***');
         end
         Inform=Inform+2;
         flag=0;
         stop = 1;
      end

      if fred0cnt > LowIts 	        % Convergence criteria 3
         if PriLev >= 1
            fprintf('*** Convergence 3, Relative function value reduction')
            fprintf(' low for %d iterations ***\n',LowIts);
         end
         flag=0;
         Inform = Inform + 4;
         stop = 1;
      end

      % Test on Lagrange multipliers

      % Dangerous to use 2nd order estimate. May ==> infinite loop
      if nract > 0 && stop && nract < n && alpha0cnt < 3 && ...
         ~isempty(H_k) && Prob.NumDiff==0
         % Check variables on bound if convergence
         % Compute search direction for free variables
         [U S V] = svd(full(H_k(idx,idx)));
         S_inv = zeros(n-nract);
         S_inv(1,1) = 1/S(1,1);
         for i = 2:n-nract
             if S(i,i) > epsRank*S(1,1)
                S_inv(i,i) = 1/S(i,i);
             end
         end
         p = V * (S_inv * (U' * (-g_k(idx))));
         % Must check how long step is possible. 
         % Otherwise 2nd order estimate might be totally wrong
         %alphaMax = 1E20;
         alphaMax = 100;
         for i=1:n-nract
             if p(i) > 0
                alphaMax = min(alphaMax,(x_U(idx(i))-x_k(idx(i)))/p(i));
             else
                if p(i) < 0
                   alphaMax = min(alphaMax,(x_L(idx(i))-x_k(idx(i)))/p(i));
                end
             end
         end
         % Compute second order estimate of Lagrange multipliers
         eta_L = lam_L(b_idx) + H_k(b_idx,idx)*alphaMax*p;
         if PriLev > 2
            xprinte(p,'p:');
            fprintf('alphaMax %16.8e\n',alphaMax)
            fprintf('2nd order Lagrange multiplier estimate I. ')
            xprinti(b_idx,'Var#:');
            xprinte(eta_L,'eta:');
         end
         % 2nd order estimate of Lagrange multiplier
         % If eta_L > 0 optimal. Use -xTol as limit
         L=double((x_U(b_idx)~=x_L(b_idx))).*eta_L < -xTol;
         if any(L)
            ix=find(L);
            stop = 0;
            Inform = 0;
            if PriLev > 2, 
               fprintf('Release variables from bound I. ');
               xprinti(b_idx(ix),'Var#:');
            end
            release=1;
            set(b_idx(ix)) = 0;
            nract = sum (set~=0);
            b_idx = find(set);
            idx = find(~set); 
         end
         lambda=eta_L;
      else
         lambda=lam_L(b_idx);
      end
      if stop && nract < n && ~isempty(H_k) && Prob.NumDiff==0 && ~issparse(H_k)
         % Check if saddle or minimum
         [V,D]=eig(H_k(idx,idx));
         V = real(V);
         D = real(D);
         [min_eig jD]=min(diag(D));
         if PriLev >=1
            fprintf('Minimal eigenvalue at stationary point %30.20f\n',min_eig);
         end
         if min_eig < -1E-12
            if eig_step || EigStepTried
               if PriLev >= 0
                  fprintf('Convergence to a saddle point. FAILURE! ');
                  fprintf('Eigenvalues:\n');
                  xprint(diag(D),'eig:');
               end
               flag=1;
               Inform=103;
            else
               stop=0;
               if PriLev >=2
                  fprintf('Saddle point. Try eigenvalue direction\n');
               end
               p=V(:,jD);
               if g_k(idx)'*p > 0 % Check descent
                  p=-p;
               end
               EigStepTried=1;
               method=method+100; % Flag that p already computed
            end
         end
      end

   %   if nract==n
   %      stop=1;
   %      flag=0;
   %      Inform=Inform+8;
   %      if PriLev >= 1
   %         disp('*** Convergence 4, Local min with all variables on bounds!');
   %         disp('*** No check has been done on the gradient.');
   %      end
   %   end
   end

   if k >= MaxIter		% Stop criteria 1
      if PriLev >= 1
         disp('*** STOP 1! Max no of iterations reached ***');
      end
      if Inform >0 && Inform < 100 
          flag=0; 
      else
          flag=1; Inform=101; 
      end
      stop = 1;
   end

   if f_k <= fLow   % Stop criteria 2
      if PriLev >= 1
         disp('*** STOP 2! Function value below Prob.f_Low   ***');
         disp('*** Restart with lower value if minimum not reached ***');
      end
      if Inform >0 && Inform < 100
          flag=0; 
      else
          flag=1; Inform=102; 
      end
      stop = 1;
   end

   if stop			% Show final results before return
      if PriLev >= 5  % REMARK! These results are displayed in the driver
         fprintf('==================================================\n');
         fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);

         xprint(x_k,'x_k:');
         xprint(set,'set:',' %2.0f',25);
         xprinte(g_k,'g_k:');
  
         if ~isempty(H_k) && Prob.NumDiff==0
            fprintf('The Hessian H_k\n');
            mPrint(H_k,'H_k:',' %10.6e',6);
         elseif ~isempty(H_k) && Prob.NumDiff ~= 0
            fprintf('Numerically computed Hessian matrix H_k\n');
            mPrint(H_k,'H_k:',' %10.6e',6);
         end

      end
      Result.Iter=k;
      Result.ExitFlag=flag;
      Result.Inform=Inform;
      Result.ExitText=ExitText(Inform);
      Result.f_k=f_k;
      Result.g_k=g_k;
      Result.B_k=B_k;
      Result.H_k=H_k;
      Result.x_k=x_k;
      Result.v_k = lambda;
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=double(x_k==x_L)+2*double(x_k==x_U);
      Result=endSolve(Prob,Result);
      return;
   end
	
release=1;
while release && ~(alpha==0 && ~isempty(xR))
%
% Determine search direction p
%
   if method == 0			% Singular Value Decomposition
      if Alg==1
         [U S V] = svd(full(H_k(idx,idx)));
      else
         [U S V] = svd(B_k(idx,idx));
      end
      S_inv = zeros(n-nract);
      if ~isempty(S)
         S_inv(1,1) = 1/S(1,1);
         if DEBUG    % HKH DEBUG output
            fprintf('METHOD 1 - Inv cond number: %30.25f\n',...
                   S(n-nract,n-nract)/S(1,1))
         end
         for i = 2:n-nract
             if S(i,i) > epsRank*S(1,1)
                S_inv(i,i) = 1/S(i,i);
             end
         end
         p = V * (S_inv * (U' * (-g_k(idx))));
      else
         p=[];
      end
   elseif method == 1		% LU-decomposition; No pivoting
      if Alg==1
         [L,U] = lu(H_k(idx,idx));
      else
         [L,U] = lu(B_k(idx,idx));
      end
      z = L\(-g_k(idx));
      p = tomsol(3, U, z);
   elseif method == 2		% LU-decomposition; with pivoting
      if Alg==1
         [L,U,P] = lu(H_k(idx,idx));
      else
         [L,U,P] = lu(B_k(idx,idx));
      end
      %z = backsubL(L, P*(-g_k(idx)));
      z = L\(P*(-g_k(idx)));
      p = tomsol(3, U, z);
   elseif method == 3		% Matlab built in inversion routine
      if Alg==1
         Hrows=size(H_k(idx,idx),1); % Add dummy column, gives QR-decomp.
         p = -([full(H_k(idx,idx)),zeros(Hrows,1)] \ g_k(idx));% Matlab does QR
      else
         Hrows=size(B_k(idx,idx),1); % Add dummy column, gives QR-decomp.
         p = -([B_k(idx,idx), zeros(Hrows,1)] \ g_k(idx)); % Matlab does QR
      end
      p = p(1:Hrows);
   elseif method == 4		% Matlab built in inversion routine
      if Alg==1
         p = -(H_k(idx,idx)\g_k(idx));
      else
         p = -(B_k(idx,idx)\g_k(idx));
      end
   elseif method == 5			% Explicit computation of inv(B_k)
      if Alg==1
         p = -(inv(H_k(idx,idx))*g_k(idx));
      else
         p = -(inv(B_k(idx,idx))*g_k(idx));
      end
   elseif method == 6			% p =-B_k*g_k, B_k = inverse of Hessian
      if Alg==1
         p = -H_k(idx,idx)*g_k(idx);
      else
         p = -B_k(idx,idx)*g_k(idx);
      end
   elseif method == 7			% p =-g_k + beta * p (Only CG)
      if cg_restart && cg_step == n+1
         cg_step=1;
         beta=0;
      end
      if k==0
         p = -g_k(idx);
      else
         p = -g_k(idx)+ beta*p_full(idx);
      end
      cg_step=cg_step+1;
   elseif method >= 100			% p already computed with eigenvectors
      method = method-100;
   end	% (if)
   release=0;
   eta_L = zeros(n,1);
   if nract > 0 && ~isempty(H_k) && ~isempty(p) && Prob.NumDiff==0
      % 2nd order estimate of Lagrange multiplier
      eta_L(b_idx) = lam_L(b_idx) + H_k(b_idx,idx)*p;
      if PriLev > 2
         fprintf('Lagrange multiplier estimate II 2nd order\n')
         xprint(eta_L(b_idx),'eta:');
      end
    
      for i=1:nract  % Check if bound constraint should be released
          j=b_idx(i);
          if lam_L(j)==0 % Only check when 1st order estimate fails
             if double((x_U(j)~=x_L(j)))*eta_L(j) < -1E-7
                set(j)=0;
                release=1;
             end
          end
      end
      if release
         nract = sum (set~=0);
         b_idx = find(set);
         idx = find(~set); % Index vector for non active variables, free vars
         if PriLev > 2
            disp('Release using 2nd order estimate')
         end
      end
   end
end
p_full=zeros(n,1);
if alpha==0 && ~isempty(xR)
   % Search in the coordinate direction of the variables to be released.
   % Last full step was a total failure
   p_full(xR)=xRp;
   p=p_full(idx);
else
   % Search step in subspace computed. Put into full space vector p_full 
   p_full(idx) = p;
end

%
% Find length of search step, alpha
%                                      1
   if isempty(idx)
      fp_0=-1;
   else
      fp_0=g_k(idx)'*p;
   end
   if k == 0
      alpha_1 = LineParam.InitStepLength;
      df=1E100;
   else
      %df = max(f_km1-f_k, 10*eps_x);
      %alpha_1 = min(1, -2*df/fp_0);
      df = max(f_km1-f_k, 10*eps_x);
      if fp_0 == 0
         alpha_1 = LineParam.InitStepLength;
      else
         alpha_1 = min(1, -2*df/fp_0);
         %alpha_1 = 1;
      end
   end
   if PriLev > 2
      fprintf('Directed derivative: %15.7e\n',fp_0);
      %fprintf('Best alpha step estimated: max(0.001,%15.7e)\n',alpha_1);
      fprintf('Best alpha step estimated: max(0.5,%15.7e)\n',alpha_1);
   end
   %alpha_1=max(1E-3,alpha_1);   % Safe guard to avoid too small step
   alpha_1=max(0.5,alpha_1);   % Safe guard to avoid too small step


   if fp_0 >= 0
      eig_step=1;
      if PriLev >= 1
         if fp_0 > 0, disp('ALARM: No descent direction.'); end
         if k > 0, fprintf('Function reduction: %15.7e  ',df);end
         fprintf('Directed derivative: %15.7e\n',fp_0);
         xprinte(p_full,'p:  ');
      end

      if ~isempty(H_k) && (~issparse(H_k) || length(idx) < 500)
         [V,D]=eig(full(H_k(idx,idx)));
         V = real(V);
         D = real(D);
         neg_eig=find(diag(D)<=0);
      else
         neg_eig=[]; % Skip using eigenvectors
      end
      if isempty(neg_eig)
      % Add negative search direction as the only possible new search direction
         if PriLev >=2
            fprintf('No negative eigenvalues. Try negative search direction\n');
         end
         f_pp = -fp_0;
         %neg_eig=1+size(V,2);
         %V=[V,-p];
         neg_eig = 1;
         V = -p;         
         D = 1;  % Dummy eigenvalue
      else
         f_pp=g_k(idx)'*V(:,neg_eig);
         iz=neg_eig(find(f_pp > 0));
         % Change sign if upward slope in eigenvector direction 
         V(:,iz)=-V(:,iz); 
         f_pp=-abs(f_pp);
         % Add negative search direction as one possible new search direction
         f_pp=[f_pp -fp_0];
         neg_eig=[neg_eig;1+size(V,2)];
         V=[V,-p];
      end

      [min_pp ix]=sort(f_pp);

      for j=1:length(f_pp);
          l=neg_eig(ix(j));
          p=V(:,l);
          alphaMax = 1E20;
          minx=alphaMax;
          xvar=0;
          for i=1:n-nract
              if p(i) > 0
                 alpha0 = (x_U(idx(i))-x_k(idx(i)))/p(i);
                 if alpha0 < alphaMax 
                    xvar=idx(i); 
                    alphaMax = alpha0;
                 end
                 %alphaMax = min(alphaMax,(x_U(idx(i))-x_k(idx(i)))/p(i));
              else
                 if p(i) < 0
                    alpha0 = (x_L(idx(i))-x_k(idx(i)))/p(i);
                    if alpha0 < alphaMax 
                       xvar=idx(i); 
                       alphaMax = alpha0;
                    end
                    %alphaMax = min(alphaMax,(x_L(idx(i))-x_k(idx(i)))/p(i));
                 end
              end
          end
          if alphaMax <= 1E-7 && j < length(f_pp)
             if PriLev >= 1
                fprintf('No move in eigenvalue direction possible!! Local minimum.\n');
                fprintf('alphaMax=%15.7e. ',alphaMax);
                fprintf('Variable %d limits the step\n',xvar);
             end
          end
          if PriLev >=3
             fprintf('Maximal step length alphaMax = %12.6f\n',alphaMax);
          end
          if alphaMax > 1E-7 % Try this eigenvalue direction
	     if PriLev >=1
                if l <= n-nract
                   fprintf('Try eigenvector %3.0f. ',l);
                   fprintf('Eigenvalue %12.5e. ',D(l,l));
                else
                   fprintf('Try negative search direction\n');
                end
                fprintf('Dir.Der.: %12.7e\n',min_pp(j));
	     end
             %g_k(idx)'*p
             p_full = zeros(n,1);
             p_full(idx) = p;
             alpha_1 = 1;
		

             LineParam.InitStepLength = alpha_1;

             LineResult = LineSearch(f, g, x_k, p_full, f_k, g_k, LineParam,...
                          alphaMax, 0, PriLev-3, Prob, varargin{:});

             alpha=LineResult.alpha;
             if plotLine
                LinePlot(f, x_k, p_full, f_k, g_k, LineParam, alphaMax,...
                         0, alpha, Prob, varargin{:});
             end
             f_k=LineResult.f_alpha;
             g_k=LineResult.g_alpha;
             alphaVec=LineResult.alphaVec;

             if PriLev > 2
                fprintf('Line search in p:')
                if alpha > 1E-5
                   fprintf(' Step length alpha = %12.6f.',alpha);
                   fprintf(' alphaMax = %12.6f\n',alphaMax);
                else
                   fprintf(' Step length alpha = %12.6e.',alpha);
                   fprintf(' alphaMax = %12.6e\n',alphaMax);
                end
             end
          else
             alpha=0;
          end
          if alpha > 1E-6 % OK
             break;
          end
      end
      f_red=0;
   else
%
%     Line search when search direction OK
%
      eig_step=0;
      if PriLev > 2
         xprinte(p_full,'p:  ');
         if PriLev > 3 || sum(set) > 0
            fprintf('Active set:');
            for j=1:n
                fprintf('%4.0f ',set(j));
            end
            fprintf('\n');
         end
         if PriLev > 3
            if alpha_1 < 0.99
               fprintf('Line search start %10.6f\n',alpha_1)
	        end
         end
      end

      alphaMax = 1E20;
      for i=1:n-nract
          if p(i) > 0
             alpha0 = (x_U(idx(i))-x_k(idx(i)))/p(i);
             if alpha0 < alphaMax 
                xvar=idx(i); 
                alphaMax = alpha0;
             end
          else
             if p(i) < 0
                alpha0 = (x_L(idx(i))-x_k(idx(i)))/p(i);
                if alpha0 < alphaMax 
                   xvar=idx(i); 
                   alphaMax = alpha0;
                end
             end
          end
      end
      f_km1 = f_k;
      if alphaMax <= 1E-14
         if PriLev >= 1
            fprintf('No move possible!! Local minimum. ');
            fprintf('alphaMax=%15.7e. ',alphaMax);
            fprintf('Variable %d limits the step\n',xvar);
         end
         alpha=0;
         f_red=0;
      else

         LineParam.InitStepLength = alpha_1;
         LineResult = LineSearch( f, g, x_k, p_full, f_k, g_k, LineParam,...
                      alphaMax, 0, PriLev-3, Prob, varargin{:});

         alpha=LineResult.alpha;
         if plotLine
            LinePlot(f, x_k, p_full, f_k, g_k, LineParam, alphaMax,...
                     0, alpha, Prob, varargin{:});
         end
         f_k=LineResult.f_alpha;
         g_k=LineResult.g_alpha;
         alphaVec=LineResult.alphaVec;

         f_red=f_km1-f_k;
         if PriLev > 2
            fprintf('Line search: alphaMax=%10.4e.',alphaMax);
            fprintf(' Func.red=%10.4e. ',f_red);
            if f_k > 0
               fprintf(' Rel.red=%10.4e. ',f_red/f_k);
            end
            if alpha > 1E-5 
               fprintf(' alpha=%12.6f. ITER %d\n',alpha,length(alphaVec));
            else
               fprintf(' Step=%12.6e. ITER %d\n',alpha,length(alphaVec));
            end
         end
      end
   end	% (if)

   if alpha < 1E-14
      alpha0cnt=alpha0cnt+1;
   else
      alpha0cnt=0;
   end

   if f_red <= eps_f * max(f_k,size_f);
      fred0cnt=fred0cnt+1;
   else
      fred0cnt=0;
   end

   if optParam.wait, pause; end;
%
% Update
%
   x_km1 = x_k;
   if CURV~=0
      xLast=x_km1;
   end
   x_k = x_k + alpha*p_full;

%  f_k and g_k already computed in LineSearch. Therefore comment next 2 lines
%  f_k = nlp_f(x_k,Prob, varargin{:});
%  g_k = nlp_g(x_k,Prob, varargin{:});

   if Alg > 1 && Alg < 6
      [B_k]=SafeUpdate(Alg,alpha,p_full,eps_x,idx,g_k,g_km1,B_k,PriLev,k);
   elseif Alg == 6 % Fletcher-Reeves CG-method
      beta=g_k(idx)'*g_k(idx) / (g_km1(idx)'*g_km1(idx));
   elseif Alg == 7 % Polak-Ribiere CG-method
      beta=(g_k(idx)-g_km1(idx))'*g_k(idx) / (g_km1(idx)'*g_km1(idx));
   elseif Alg == 8 % Fletcher conjugate descent CG-method
      beta=-(g_k(idx)'*g_k(idx)) / (g_km1(idx)'*p_full(idx));
      %fprintf('beta %20.10f\n',beta);
   end

   g_km1 = g_k;
   k = k+1;
end	% (while)

% =========================================
function [SolvAlg, Alg, method] = NameAlg(Alg, method)
% =========================================
nargin;

switch Alg
case 1
   switch method
   case 0
       SolvAlg='Newtons method using SVD';
   case 1
       SolvAlg='Newtons method using LU';
   case 2
       SolvAlg='Newtons method using LU with pivoting';
   case 3
       SolvAlg='Newtons method using built in QR';
   case 4
       SolvAlg='Newtons method using built in inversion';
   case 5
       SolvAlg='Newtons method using explicit inverse';
   otherwise
       SolvAlg='Newtons method using SVD';
       method=0;
   end
case 2
   SolvAlg='Safeguarded BFGS';
case 3
   SolvAlg='Safeguarded BFGS (Hessian update)';
case 4
   SolvAlg='Safeguarded DFP';
case 5
   SolvAlg='Safeguarded DFP (Hessian update)';
case 6
   if method == 1
      SolvAlg='Fletcher-Reeves CG with reset';
   else
      SolvAlg='Fletcher-Reeves CG';
   end
case 7
   if method == 1
      SolvAlg='Polak-Ribieres CG with reset';
   else
      SolvAlg='Polak-Ribieres CG';
   end
case 8
   if method == 1
      SolvAlg='Fletcher conjugate descent CG-method with reset';
   else
      SolvAlg='Fletcher conjugate descent CG-method';
   end
otherwise
   Alg=2;
   SolvAlg='Safeguarded BFGS';
end

% =========================================
function [B_k]=SafeUpdate(Alg,alpha,p_full,eps_x,idx,g_k,g_km1,B_k,PriLev,k)
% =========================================
switch Alg
   case 2 % Safeguarded BFGS update of the inverse Hessian
      z=alpha*p_full(idx); % Check if we have converged (or short alpha step)
      if norm(z) > eps_x
         y=g_k(idx)-g_km1(idx);
         By=B_k(idx,idx)*y;
         zy=z'*y;
         if PriLev > 2
            fprintf('BFGS Hessian update iter %3.0f\n',k);
            fprintf('zy %10.6e\n',zy)
         end
         if zy < 1E-13 
            if PriLev > 2
               fprintf('BFGS Inverse Hessian update step dangerous!\n');
               fprintf('max(z)/zy %20.10e\n',max(z)/zy);
               fprintf('min(z)/zy %20.10e\n',min(z)/zy);
            end
         end
         if zy ~= 0  
            B_k(idx,idx)=B_k(idx,idx)+((1+(y/zy)'*By)/zy)*z*z'-(z*By'+By*z')/zy;
         end
         if PriLev > 3
            fprintf('B_k matrix after BFGS update\n');
            PrintMatrix(B_k,'B_k:')
            fprintf('B_k matrix eigenvalues:\n');
            xprint(eig(B_k),'eig:');
         end
      end
   case 3 % Safeguarded BFGS update of approximate Hessian
      z=alpha*p_full(idx); % Check if we have converged (or short alpha step)
      if norm(z) > eps_x
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
            if PriLev > 2
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
   case 4 % Safeguarded DFP update of approximate inverse Hessian
      z=alpha*p_full(idx); % Check if we have converged (or short alpha step)
      if norm(z) > eps_x
         y=g_k(idx)-g_km1(idx);
         By=B_k(idx,idx)*y;
         yBy=y'*By;
         zy=z'*y;
         if zy < 0.2 * yBy
            theta=0.8*yBy / (yBy - zy);
            w=theta * z + (1-theta)*By;
            zw=y'*w;
         else 
            w=z;
            zw=zy;
         end
         if PriLev > 2
            fprintf('DFP Inverse Hessian update iter %3.0f\n',k);
            fprintf('zw %10.6e yBy %10.6e \n',zw,yBy)
         end
         if zw < 1E-13 | yBy < 1E-13
            if PriLev > 2
               fprintf('DFP Inverse Hessian update step dangerous!\n');
               fprintf('max(w)/zw %20.10e\n',max(w)/zw);
               fprintf('min(w)/zw %20.10e\n',min(w)/zw);
               fprintf('max(By)/yBy %20.10e\n',max(By)/yBy);
               fprintf('min(By)/yBy %20.10e\n',min(By)/yBy);
            end
         end
         if zw == 0  
            if yBy ~=0
               B_k(idx,idx)=B_k(idx,idx)-(By/yBy)*By';
            end
         elseif  yBy == 0
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w';
         else
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w'-(By/yBy)*By';
         end
         if PriLev > 3
            fprintf('B_k matrix after Inverse DFP update\n');
            PrintMatrix(B_k,'B_k:')
            fprintf('B_k matrix eigenvalues:\n');
            xprint(eig(B_k),'eig:');
         end
      end
   case 5 % Safeguarded DFP update of the approximate Hessian
      z=alpha*p_full(idx); % Check if we have converged (or short alpha step)
      if norm(z) > eps_x
         y=g_k(idx)-g_km1(idx);
         Bz=B_k(idx,idx)*z;
         zy=z'*y;
         if PriLev > 2
            fprintf('DFP Hessian update iter %3.0f\n',k);
            fprintf('zy %10.6e\n',zy)
         end
         if zy < 1E-13 
            if PriLev > 2
               fprintf('DFP Hessian update step dangerous!\n');
               fprintf('max(y)/zy %20.10e\n',max(y)/zy);
               fprintf('min(y)/zy %20.10e\n',min(y)/zy);
            end
         end
         if zy ~= 0  
            B_k(idx,idx)=B_k(idx,idx)+((1+(z/zy)'*Bz)/zy)*y*y'-(y*Bz'+Bz*y')/zy;
         end
         if PriLev > 3
            fprintf('B_k matrix after DFP update\n');
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
     Text = 'Relative f(x) reduction low for 10 iterations';
   case 5
     Text = 'Iterations close, relative f(x) reduction low';
   case 6
     Text = 'Proj grad small, f(x) reduction low for 10 iterations';
   case 7
     Text = str2mat('Iterations close, projected gradient small' ...
                   ,'f(x) reduction low for 10 iterations');
   case 101
     Text = 'Maximal number of iterations reached';
   case 102
     Text = str2mat('Function value below AbsF' ...
                   ,'Restart with lower value if minimum not reached');
   case 103
     Text = 'Convergence to saddle point (eigenvalues computed)';
end
if Inform < 100
   Text=str2mat('Optimal solution found',Text);
end

% MODIFICATION LOG:
%
% 980618  mbk   Second criteria for releasing variables aborted. Now, variable
%               can't be released if activated in previous iteration unless
%               alpha=0.
%               Parameter not_releaseall not used, deleted.
%               Set V=[-p] instead of V=[V,-p] when isempty(H).
%               Parameter testconv not used, deleted.
%               Change of Alg=1 and Alg=4 since mixed up in approximate Hessian
%               update.
%
% 980622  mbk   Set Result.v_k to lam_L or eta_L.
%               Set p equal to the eigenvector corresponding to the most
%               negative eigenvalue, p=V(:,j). Line ~= 411.
%
% 980624  mbk   Double computation of 'df' and 'alpha_1' aborted with '%'
%               operator. Line ~= 581
%
% 980918  hkh   Changed method to range [0,8] instead of [1,9]
%               Use of structure optParam, instead of optPar vector
% 980919  hkh   Call new LineSearch instead of linesrch, change variable names
% 980920  hkh   Bug, D=1 must be set if no eigenvalues possible.
% 981005  mbk   Changes in comments concerning type of convergence.
% 981005  hkh   Change bTol to xTol. Better convergence tests.
%               New call to LineSearch. Changing lineplot to LinePlot.
%               New flags linePlot and DEBUG. Changed inform to Inform
% 981011  hkh   Added Prob to Result struct as default
% 981013  hkh   Added call to iniSolve and endSolve
% 981017  hkh   Deleted the gradient as argument to LinePlot
% 981026  hkh   Solver is name of solver. SolverAlgorithm is description
%               Avoid infinite loop in eigenvalue directions, EigStepTried
% 981026  hkh   Put Prob.f_Low in optParam.fLow to LineSearch
% 981027  hkh   Use Prob.Solver.Alg to define Alg
% 981028  hkh   Put f_0 = f(x_0) as field in Result
% 990306  hkh   Test on FUNCS.H instead of p_H, to see if Hessian defined
% 990308  hkh   Add comments about the setting of Prob.Solver.Alg
%               Safeguard against empty H_k, switch Newton to BFGS in such case
% 990311  hkh   Safeguard against empty initial g_k, use numerical differences.
% 000904  hkh   Change default behaviour to either Newton or BFGS
% 000917  hkh   Add ExitText
% 000919  hkh   Check on Prob.FUNCS.g if user has not given derivative
% 000923  hkh   Revise parameter handling
% 001106  hkh   Use method from Prob.Solver.Method
% 001206  hkh   Take real of eigval/vect, safe against numerical fuzz
% 010223  hkh   Change alphaMax to 100 instead of 1E20, more safe
% 020409  hkh   Use Prob.nState and Prob.Mode
% 020819  ago   Modify for Matlab 6.5 logical handling
% 031201  hkh   Change comments for AD, new fields ADObj
% 040111  hkh   Change call to inisolve
% 040916  med   logical mtimes fixed
% 040920  frhe  Logicals casted to doubles to prevent mult's with logicals.
% 041018  hkh   Allow LineParam.sigma [], set 0.9 if so
% 060712  med   Updated help
% 060814  med   FUNCS used for callbacks instead
% 060818  hkh   Use Prob.f_Low instead of optParam.eps_absf for target f test
% 061212  med   ADMAT removed
% 070728  hkh   Change print level for "No move" to >= 1 instead of > -1
% 070728  hkh   Avoid "No move" output for negative search direction, change text
% 080604  hkh   Safeguard check on NumDiff, if []
% 080606  hkh   % SICK Prob.ConsDiff = 0, avoid output in case misuse with c(x)
% 090720  med   Help updated, code cleaned