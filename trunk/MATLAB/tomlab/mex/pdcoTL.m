% TOMLAB PDCO Solver, calls pdco.m:
%
% Primal-Dual Barrier Method for Convex Objectives
%
% Developed by Stanford Systems Laboratory (SOL)
%
% function Result = pdcoTL(Prob)
%
% Tomlab pdco solves linearly constrained problems with convex objectives:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds
%    b_L <= A*x  <= b_U, linear constraints
%
% The linear constraints are transformed to the SOL PDCO problem form:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds
%           A*x   = b  , linear constraints
%
% pdco actually solves the following problem:
%
%    minimize    f(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
%
% where
%    D1, D2 are positive-definite diagonal matrices defined from d1, d2.
%           In particular, d2 indicates the accuracy required for
%           satisfying each row of Ax = b.
%           See the help for pdco.m for a detailed discussion of D1 and D2
%           In pdco.m the objective is called phi(x), not f(x)
%           and bl == x_L, bu == x_U
%
% ---------------------------------------------------------------------
%               ---  PDCO parameters ---
%
%               PDCO takes a number of input parameters as defined in the
%               file pdcoSet.m (see help for pdcoSet, and input help below).
%               The parameters are fields in the structure options
%
%               This structure is possible to send as Prob.SOL.pdco
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%  x_0          Initial x vector, used if non-empty.
%
%  x_L, x_U     Bounds on variables. x_L(k) = x_U(k) if x(k) fixed.
%  b_L, b_U     Bounds on linear constraints.
%               For equality constraints, set b_L(k) == b_U(k) if k equality.
%
%  A            Matrix of coefficients for the linear constraints.
%
%  PriLevOpt    Print level in pdco solver. If > 0 prints summary information.
%
%  SOL          Structure with special fields for SOL solvers.
%  ===          The following fields are used:
%
%    d1         D1=diag(d1), if d1 n-dimensional
%               D1=diag(d1*ones(n,1)), if d1 is a scalar (default 1E-4)
%    d2         D2=diag(d2), if d2 m-dimensional
%               D2=diag(d2*ones(m,1)), if d2 is a scalar (default 1E-4)
%    y0         Initial dual parameters for linear constraints (default 0)
%    z0         Initial dual parameters for simple bounds (default 1/N)
%    xsize      Estimate of the biggest x at the solution. (default 1)
%    zsize      Estimate of the biggest z at the solution. (default 1)
%               xsize,zsize are used to scale (x,y,z).  Good estimates
%               should improve the performance of the barrier method
%    pdco       Structure with the same fields as pdcoSet defines
%               Fields in Prob.SOL.pdco
%      MaxIter     Maximum iterations of the primal-dual barrier method.
%                  Most problems should solve within 30 PDitns. Default 500
%      FeaTol      Accuracy for satisfying Ax + D2 r = b, A'y + z = gobj
%                  and x - x1 = bl, x + x2 = bu, where x1, x2 > 0.
%                  1e-6 is typically small enough.
%                  1e-5 may be acceptable also.
%      OptTol      Accuracy for satisfying x1.*z1 = 0, x2.*z2 = 0,
%                  where z = z1 - z2 and z1, z2 > 0.
%                  Typically the same as Featol.
%      StepTol     (between 0 and 1): Controls how close each an
%                  x or z may be to reaching a bound at each step.
%                  For safety, should not be bigger than 0.99 (say)
%                  for nonlinear problems.
%      StepSame    1 (true) if stepx and stepz should be the same
%                  (gives true Newton method for nonlinear problems);
%                  0 (false) if stepx and stepz may be different
%                  (gives better performance for linear programs).
%      x0min       Min distance between x0 and bl or bu  AFTER SCALING.
%                  1.0 is about right for cold starts.
%                  0.1 or 0.01 may be ok if x0, z0 are "good".
%      z0min       Min distance between abs(z0) and zero AFTER SCALING,
%                  when z0 is used to initialize z1 > 0 and z2 > 0.
%                  Typically the same as x0min.
%      mu0         Initial mu (ABSOLUTE VALUE) for solving scaled problem.
%      Method      Specifies how each search direction (dx,dy,dz1,dz2)
%                  should be computed.  Several methods exist for
%                  experimental purposes.  If A has fewer rows than columns
%                  (m < n) we usually solve for dy first (most of the work)
%                  and then get the other components cheaply.
%
%  Method  Solve for  Using                                         Implemented?
%     1       dy      Sparse Cholesky on (A D^2 A' + D2^2 I).           Yes
%     2       dy      Sparse QR on corresponding least-squares problem  Yes
%     3       dy      LSQR on least-squares problem                     Yes
%
%    11       dx      Sparse Cholesky on (D A'A D  + D2^2 I).           No
%    12       dx      Sparse QR on corresponding least-squares problem  No
%    13       dx      LSQR on least-squares problem                     No
%
%    21    dx,dy      Sparse LU on 2x2 KKT-type system                  No
%    23    dx,dy      SYMMLQ on same system                             No
%
%    31    dx,dy,dz1  Sparse LU on 3x3 system (only for problems with   No
%                     with vanilla bounds: 0 < x < inf).
%                     This is a HUGE system, but it is relatively
%                     well-conditioned compared to the other systems.
%
%    41 dx,dy,dz1,dz2 Sparse LU on 4x4 system with general bounds.      Yes
%                     This is an even HUGER system.
%
%    4 = Tomlab Tlsqr, special PDCO interface avoiding any callbacks    Yes
%    5 = Tomlab Tlsqr, special interface, with callbacks                Yes
%    6 = Tomlab Tlsqr, atol algorithm in Matlab, with callbacks         Yes
%
% Option 5 and 6 are slower than 4, and normally not to be used
%
%    If A is an explicit sparse matrix, all methods are applicable.
%     1 is usually best (e.g. for LPs).
%     2 may be more reliable; it's there for checking.
%     3 is sometimes more efficient (e.g. for entropy problems).
%       Diagonal preconditioning is possible.
%
%    If A is an operator,
%     3 must be used.  Diagonal preconditioning is not possible.
%
%    Notes:
%    Method =  1      On ENTROPY.big, symamd can find an ordering
%                     but chol never finishes.
%    Method = 41      On ENTROPY.big, symamd never finishes.
%
% The following options control LSQR when Method = 3,4,5,6:
%
% LSQRMaxIter * min(m,n) is the maximum LSQR (CG) iterations.
% LSQRatol1   is the starting value of the LSQR accuracy
%                     tolerance "atol" (if LSmethod = 3).
%                     1e-3 or 1e-4 sometimes works.
%                     1e-8 may be needed for LPs.
%                     In general, if max(Pinf,Dinf,Cinf) doesn't decrease
%                     every iteration, set atol1 to a smaller value.
% LSQRatol2   is the smallest value atol is reduced to.
% LSQRconlim  shuts LSQR down early if its matrix is ill-conditioned.
%
% wait = 0    means pdco should proceed to solve the problem. (DISABLED)
%              = 1    means pdco should pause to allow interactive resetting
%                     of some of the parameters.
%
% ------------------------------------------------------------------------------
%
% OUTPUT:
% Result        Structure with optimization results
%
%   f_0         Function value at start, x = x_0.
%   f_k         Function value at optimum.
%   g_k         Gradient of the function at the solution.
%   H_k         Hessian of the function at the solution, diagonal only.
%
%   x_k         Solution vector.
%   x_0         Initial solution vector.
%
%   xState      State of variables. Free == 0; On lower == 1; On upper == 2;
%               Fixed == 3;
%   bState      State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   v_k         Lagrangian multipliers (orignal bounds + constraints).
%
%   y_k         Lagrangian multipliers (for bounds + dual solution vector)
%               The full [z;y] vector as returned from pdco, including slacks
%               after rewriting constraints:
%               -inf < b_L < A*x < b_U < inf;  to Ax + s = b, s slack variables
%               For non-inf lower AND upper bounds the corresponding slack variable
%               tells which of the bound (if any) that are binding.
%
%   ExitFlag    Tomlab Exit status from pdco MEX.
%   Inform      pdco information parameter: 0 = Solution found;
%               1 = Too many iterations; 2 = Linesearch failed too often.
%
%   Iter        Number of iterations.
%   FuncEv      Number of function evaluations.
%   GradEv      Number of gradient evaluations.
%   HessEv      Number of Hessian  evaluations.
%
%   Solver           Name of the solver (pdco).
%   SolverAlgorithm  Description of the solver.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Jan 20, 2003.   Last modified May 26, 2007.

function Result=pdcoTL(Prob)

%#function pd_fgh

if nargin < 1, error('pdcoTL needs the Prob structure as input');return;end

Prob.solvType = 3;

Prob = iniSolve(Prob,3,2,0);

probType = Prob.probType;

Result=ResultDef(Prob);
Result.Solver = 'pdco';

%LargeScale = DefPar(Prob,'LargeScale',1);

[bl,bu,n,nnLin,nnCon] = defblbu(Prob,inf,1);

if nnCon > 0
   fprintf('ERROR - PDCO only solves linearly constrained problems\n');
   Result.ExitText = 'PDCO does not allow nonlinear constraints\n';
   Result.ExitFlag = 10;
   Result=endSolve(Prob,Result);
   return
end

LP = checkType('lp');
LLS = checkType('lls'); % Linear LS
NLLS = checkType('ls'); % LS with only bounds

PriLev=DefPar(Prob,'PriLevOpt',0);
%HKH set 1 now for debugging
PriLev = 0;
Prob.PriLevOpt = PriLev;

%if isempty(Prob.HessPattern)
%   % Only the diagonal of the Hessian is needed
%   Prob.HessPattern = speye(n);
%end
%HKH Always set to Hessian pattern to the diagonal, pdco only uses the diagonal 

Prob.HessPattern = speye(n);

options.MaxIter      =  1000;
options.FeaTol       =  1e-6;
options.OptTol       =  1e-6;
if any(probType == [LP,LLS])
   options.StepSame     =       0;  % 1 for stepx == stepz (NLPs)
   options.StepTol      =  0.9999;
   options.LSQRMaxIter  =    10.0;
   options.LSQRatol1    =    1e-8;
   options.LSQRatol2    =   1e-15;
else
   options.StepSame     =     1;  % 1 for stepx == stepz (NLPs)
   options.StepTol      =  0.99;
   options.LSQRMaxIter  =  10.0;
   options.LSQRatol1    =  1e-6;
   options.LSQRatol2    = 1e-13;  % 
end
options.LSQRconlim   = 1e+12;  % Somewhere between e+8 and e+16
if Prob.WarmStart | ~isempty(Prob.x_0)
   options.x0min        =   0.1;  % 1.0 | 0.1 for cold | warm starts?
   options.z0min        =   0.1;  % 1.0 | 0.1 for cold | warm starts?
else
   options.x0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
   options.z0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
end
options.mu0          = 1e-1;  % 1.0 starts near central path
% Default for non-Tomlab use
% defoptions.Method       =     3;  % 3 = computed dy using LSQR
% Default for use in Tomlab
options.Method       =   4;  % 4 = computed dy using Tlsqr
options.wait         =     0;
options.NOTE         = 'LSQRMaxIter is scaled by the matrix dimension';

if isfield(Prob.SOL,'pdco')
   UserOpt             = Prob.SOL.pdco;
   options.MaxIter     = DefPar(UserOpt,'MaxIter',options.MaxIter);
   options.FeaTol      = DefPar(UserOpt,'FeaTol',options.FeaTol);
   options.OptTol      = DefPar(UserOpt,'OptTol',options.OptTol);
   options.StepTol     = DefPar(UserOpt,'StepTol',options.StepTol);
   options.StepSame    = DefPar(UserOpt,'StepSame',options.StepSame);
   options.x0min       = DefPar(UserOpt,'x0min',options.x0min);
   options.z0min       = DefPar(UserOpt,'z0min',options.z0min);
   options.mu0         = DefPar(UserOpt,'mu0',options.mu0);
   options.Method      = DefPar(UserOpt,'Method',options.Method);
   options.LSQRMaxIter = DefPar(UserOpt,'LSQRMaxIter',options.LSQRMaxIter);
   options.LSQRatol1   = DefPar(UserOpt,'LSQRatol1',options.LSQRatol1);
   options.LSQRatol2   = DefPar(UserOpt,'LSQRatol2',options.LSQRatol2);
   options.LSQRconlim  = DefPar(UserOpt,'LSQRconlim',options.LSQRconlim);
   % options.wait = DefPar(UserOpt,'wait',options.wait);
   % options.NOTE = DefPar(UserOpt,'NOTE',options.NOTE);
end

D1 = DefPar(Prob.SOL,'d1',1E-4);
D2 = DefPar(Prob.SOL,'d2',1E-4);


x_L = bl(1:n);
x_U = bu(1:n);

if nnLin > 0
   b_L = bl(n+1:n+nnLin);
   b_U = bu(n+1:n+nnLin);
   b_D = b_U-b_L;
   if all(b_D == 0)
      % Pure equality problem
      A = Prob.A;
      b = b_U;
   else
      % Rewrite to equality problem adding slacks
      b   = b_U;
      s_L = zeros(nnLin,1); 
      s_U = zeros(nnLin,1); 
      s_U(isinf(b_L)) = inf; 
      ixU = find(isinf(b_U));
      s_L(ixU) = -inf; 
      b(ixU)   = b_L(ixU);
      ixI = find(~isinf(b_D));
      s_U(ixI) = b_D(ixI);

      % New problem:
      A = [Prob.A, speye(nnLin,nnLin)];
      x_L =  [x_L; s_L];
      x_U =  [x_U; s_U];
   end
else
   % Create a dummy constraint if no linear constraints present
   %
   % x(ix) + x(n+1) = x_U(ix), where ix is first non-inf upper bound
   % if all upper bounds are inf
   % x(ix) - x(n+1) = x_L(ix), where ix is first non-inf lower bound

   ix = find(~isinf(x_U));
   if ~isempty(ix)
      ix        = ix(1);
      A         = zeros(1,n+1);
      A(1,ix)   = 1;
      A(1,n+1)  = 1;
      b_L       = -inf;
      b_U       = x_U(ix);
      b         = b_U;
   else
      ix = find(~isinf(x_L));
      ix        = ix(1);
      A         = zeros(1,n+1);
      A(1,ix)   = 1;
      A(1,n+1)  = -1;
      b_L       = x_L(ix);
      b_U       = inf;
      b         = b_L;
   end
   x_L =  [x_L; 0];
   x_U =  [x_U; inf];
end

N = length(x_L);

% Starting point
if isempty(Prob.x_0)
   %x_0 = ones(N,1)/N;   
   if nnLin > 0 & length(x_L) > n
      % Add initial slack values
      x_0 = [ones(n,1)/N;b-A(:,1:n)*ones(n,1)];   
   elseif nnLin > 0 & length(x_L) == n
      x_0 = ones(n,1)/N;
   elseif nnLin == 0
      % Add initial value for the dummy constraint slack
      x_0 = [ones(n,1)/N;b-1/n];   
   end
else
   x_0  = DefPar(Prob,'x_0',ones(n,1)/N);
   if nnLin > 0 & length(x_L) > n
      % Add initial slack values
      x_0  = [x_0;b-A(:,1:n)*x_0];
   elseif nnLin == 0
      % Add initial value for the dummy constraint slack
      x_0  = [x_0;b-A(:,1:n)*x_0];
   end
   % if isempty(xsize), xsize = max(max(abs(x_0)),1); end
end
xsize = DefPar(Prob.SOL,'xsize',[]);
%if isempty(xsize), xsize = 1/N; end
% Mike suggests default xsize, zsize to be 1, 1/N good for entropy only
if isempty(xsize), xsize = 1; end
if xsize == 0, xsize = 1/N; end

% [x_L x_0 x_U]
% Safe-guard starting point
x_0 = max(x_L(1:N),min(x_0,x_U(1:N)));

y0 = DefPar(Prob.SOL,'y0',[]);

if isempty(y0) 
   y0 = zeros(length(b),1);
   % y0 = (b-A*x_0);
end
if length(y0) < nnLin
   y0 = [y0;zeros(nnLin-length(y0),1)];
end
% y0 = y0./(D2.^2);


z0 = DefPar(Prob.SOL,'z0',[]);
if isempty(z0)
   % z0 = ones(N,1)/N;   
   z0 = zeros(N,1);   
   z0(x_0==x_U) = -1/N;
   z0(x_0==x_L) = 1/N;
end
if length(z0) < N
   z0 = [z0;ones(N-length(z0),1)];
end
zsize = DefPar(Prob.SOL,'zsize',[]);
%if isempty(zsize), zsize = 1/N; end
if isempty(zsize), zsize = xsize; end

% Mike suggests default xsize, zsize to be 1, 1/N good for entropy only
% if isempty(zsize), zsize = 1; end
% if isempty(zsize), zsize = max(max(z0),1); end
if zsize == 0, zsize = 1; end


%ix = find(x_L~=x_0 & x_U~=x_0);
%options.x0min = min(min(abs((x_L(ix)-x_0(ix))/xsize),abs((x_U(ix)-x_0(ix))/xsize)))
%options.z0min = options.x0min/10;
%x0min=options.x0min
%z0min=options.z0min

Result.f_0  = nlp_f(x_0(1:n),Prob);

switch options.Method
   case 1
      Result.SolverAlgorithm = 'TOMLAB PDCO, using Cholesky';
   case 2
      Result.SolverAlgorithm = 'TOMLAB PDCO, using QR';
   case 3
      Result.SolverAlgorithm = 'TOMLAB PDCO, using SOL Matlab LSQR';
   case 41
      Result.SolverAlgorithm = 'TOMLAB PDCO, using sparse LU';
   case 5
      Result.SolverAlgorithm = 'TOMLAB PDCO, Tlsqr with callbacks';
   case 6
      Result.SolverAlgorithm = 'TOMLAB PDCO, Iterative Tlsqr with callbacks';
   otherwise  % Case 4 default
      Result.SolverAlgorithm = 'TOMLAB PDCO, using Tomlab Tlsqr Special Mex';
end

%if 0
%ProbS = Prob;
%ProbS.A    = A;
%ProbS.b_L  = b;
%ProbS.b_U  = b;
%ProbS.N    = N;
%ProbS.x_0  = x_0;
%ProbS.x_L  = x_L;
%ProbS.x_U  = x_U;
%ProbS.FUNCS.f = 'ls_ff';
%ProbS.FUNCS.g = 'ls_gg';
%ProbS.FUNCS.H = '';
%rr=tomRun('snopt',ProbS,2);
%end

%y0 = zeros(nnLin,1);
%z0 = ones(N,1)/N;
%z0 = x_0;
%z0 = ones(N,1);
%options.mu0          =  0.01;  % < 1.0 better than 1.0?
%options.mu0          =  1E-5;  % < 1.0 better than 1.0?
%options.x0min       = 0.1;
%options.z0min       = 0.1;
%xsize = 1.0;
%zsize = 1.0;
%options.FeaTol       =  1e-6;
%options.OptTol       =  1e-6;
%D1 = 1E-4;
%D2 = 1E-4;
%D2 = [ones(n,1); 1E-3*ones(nnLin,1)];
%options.LSQRatol1    = 1e-4;
%options.LSQRatol2    = 1e-8;   % Let LSQR solve loosely to start with
%options.LSQRMaxIter  = 30;
%options.StepSame     =     0;  % 1 for stepx == stepz (NLPs)
if PriLev > 0
   options
end

Inform =2;
Iter = 0;
i    = 0;

r = b-A*x_0;
rnorm = norm(r);
if PriLev > 0
   rnorm
end
%alfa = 0.1;
%alfa = N;
alfa = 1;
feasible = 0;
if probType == LP
   iMax = 2;
else
   iMax = 6;
end

while Inform == 2 & Iter < options.MaxIter & i < iMax
   [x,y,z,Inform,PDitns,CGitns,time] = ...
      pdco('pd_fgh',A,b,x_L,x_U,D1,D2,options,x_0,y0,z0,xsize,zsize, Prob );
   Iter = Iter + PDitns;
   x_0 = x;
   %xsize = max(x_0);
   %zsize = xsize;
   %y0  = y;
   % z0  = z;
   y0  = zeros(nnLin,1);
   alfa = xsize;
   %z0  = ones(N,1);
   %z0  = ones(N,1)/N;
   z0  = alfa*ones(N,1)/N;
   %z0 = zeros(N,1);   
   %z0(x_0==x_U) = -alfa;
   %z0(x_0==x_L) = alfa;
   %z0 = 2*alfa*rand(N,1)-alfa;
   %z0 = 2*alfa*rand(N,1);
   i   = i+1;
   r = b-A*x_0;
   rnorm = norm(r);
   if PriLev > 0
      rnorm
   end
   if rnorm < 1E-7 & ~feasible
      options.mu0 = 0.001;
      feasible = 1;
   elseif feasible
      if PDitns < 5
         options.mu0 = min(1,options.mu0*10); % 1.0 starts near central path
      else
         %options.mu0 = min(1,options.mu0*2); % 1.0 starts near central path
         options.mu0 = min(1,options.mu0/2); % 1.0 starts near central path
      end
   elseif probType == LP
      % LP, not feasible
      if i==1
         options.mu0 = 0.1;
      else
         options.mu0 = min(1,options.mu0*2); % 1.0 starts near central path
      end
      %options.StepSame     =     1-options.StepSame; 
      options.StepSame     =     1;
      %options.mu0 = 0.9;
      %options.mu0 = 0.001;
      %options.mu0 = 0.5;
      %xsize  = 1/N;
      %ysize  = 1/N;
      %options.LSQRatol1 = 1e-6;  % Let LSQR solve loosely to start with
      %z0  = rand(N,1)/N;
      %z0  = zeros(N,1);
   else
      % non-LP, not feasible
      options.mu0 = min(1,options.mu0*2); % 1.0 starts near central path
   end
   if feasible
      y0  = y;
   end
   %if i == 2
      %z0 = zeros(N,1);   
      %y0  = y;
   %end
end

%plot(b-A*x)
%sum(max(0,b-A*x))

Result.Inform    = Inform;
Result.Iter      = PDitns;
Result.MinorIter = CGitns;

Result.x_k  = x(1:n);
Result.f_k  = nlp_f(x(1:n),Prob);
Result.x_0  = x_0(1:n);
Result.y_k  = [z;y];

if ~isempty(Prob.FUNCS.g),Result.g_k  = nlp_g(x(1:n),Prob);end
if ~isempty(Prob.FUNCS.H),Result.H_k  = nlp_H(x(1:n),Prob);end

if nnLin > 0
   Result.v_k  = [z(1:n);y];
else
   Result.v_k  = z(1:n);
end

global n_f n_g n_H
Result.FuncEv   = n_f;
Result.GradEv   = n_g;
Result.HessEv   = n_H;

optParam = Prob.optParam;

if(nnLin>0)
  %Ax = Prob.A*x_k;
  Result.Ax = Prob.A * x(1:n);
else
  Result.Ax = [];
end

% Must check x and Ax with looser tolerance ???
%eps_x = DefPar(options,'OptTol',1.73E-6);
%bTol = DefPar(options,'LSQRatol',1E-8);

if nnLin > 0
   Result = StateDef(Result, x(1:n), Result.Ax, [], ...
                     options.OptTol, options.FeaTol, [], bl, bu,1);
else
   Result = StateDef(Result, x(1:n), [], [], ...
                     options.OptTol, options.FeaTol, [], bl, bu,1);
end

switch Result.Inform
   case 1
     ExitFlag=1;   % Too many iterations
   %case 1
   %  ExitFlag=2;  % Unbounded
   %case {2,3,4}
   %  ExitFlag=4;  % Infeasible
   case {2}
     ExitFlag=3;   % Rank problem
   %case {9,10,7}
   %  ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end
Result.ExitFlag = ExitFlag;

switch(Result.Inform)
case 0,
   Result.ExitText = 'Solution found';
case 1,
   Result.ExitText = 'Too many iterations';
case 2,
   Result.ExitText = 'Linesearch failed too often';
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 030120 hkh  Written
% 030126 hkh  Input and output completed
% 030406 hkh  b=[] needed if no linear constraints
% 030818 ango Fixed bug when Prob.A is one row
% 031121 hkh  Removed LSproblem, and changed LSmethod to Method
% 031204 ango Changed pd_fgH reference to pd_fgh
% 031208 hkh  Mike S suggested default xsize,zsize be 1 instead of 1/N
% 040109 hkh  Revised for v4.2, iniSolve, endSolve, no ProbCheck
% 040803 med  Pragmas added for MATLAB Compiler
% 041201 hkh  Call endSolve before error return
% 041202 hkh  Revise calls to defblbu and StateDef, use output from defblbu
% 041222 med  Safeguard added for x_0
% 050109 med  Safeguard corrected
% 050726 med  x_0 safeguard corrected again
% 060801 hkh  Revise b_L<=Ax<=b_U to Ax=b reformulation, avoid change of A
% 060801 hkh  Revise parameter handling
% 060802 hkh  Implement adaptive strategy to handle linesearch failure (common)
% 060814 med  FUNCS used for callbacks instead
% 070526 med  Output removed