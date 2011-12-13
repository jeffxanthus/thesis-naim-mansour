% TOMLAB PDSCO Solver, calls pdsco.m: 
%
% Primal-Dual Barrier Method for Separable Convex Objectives
%
% Developed by Stanford Systems Laboratory (SOL)
%
% function Result = pdscoTL(Prob)
%
% TOMLAB pdsco solves linearly constrained problems with convex objectives:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds, with x nonnegative
%    b_L <= A*x  <= b_U, linear constraints
%
%    NOTE!!! Only x_L == 0, and x_U == inf is allowed
%
% The linear constraints are transformed to the SOL PDSCO problem form:
%
%    min f(x)
%
%    subject to
%
%       x   >= 0, variable bounds
%       A*x  = b, linear constraints
%
% pdsco actually solves the following problem:
%
%    minimize    f(x)  +  1/2 norm(gamma*x)^2  +  1/2 norm(r/delta)^2
%      x,r
%    subject to  Ax + r = b,    x > 0,    r unconstrained.
%
% where
% f(x)   is a smooth separable convex function (possibly linear);
% gamma  is a primal regularization parameter (typically small but may be 0);
% delta  is a dual regularization parameter (typically small or 1; must be >0);
% With positive gamma and delta, the primal-dual solution (x,y,z) is 
% bounded and unique.
%
%           See help pdsco.m for a detailed discussion of gamma and delta.
%           In pdsco.m the objective is called phi(x), not f(x).
%
% ---------------------------------------------------------------------
%               ---  PDSCO parameters ---
%
%               PDSCO takes a number of input parameters as defined in the 
%               file pdscoSet.m (see help).
%               The parameters are fields in the structure options.
%
%               This structure is possible to send as Prob.SOL.pdco.
%
%               Some of the fields in Prob.SOL.pdco also have counterparts
%               as standard parameters in Tomlab. If these parameters are
%               changed from their default values (as defined by pdscoSet),
%               then they are used, e.g. 
%               Prob.optParam.MaxIter (=Prob.SOL.pdco.MaxIter).
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
%  PriLevOpt    Print level in pdsco solver. If > 0 prints summary information.
%
%  SOL          Structure with special fields for SOL solvers.
%  ===          The following fields are used:
%
%    pdco       Structure with exactly the same fields as pdscoSet defines
%    y0         Initial dual parameters for linear constraints (default 0)
%    z0         Initial dual parameters for simple bounds (default 1/N)
%    xsize      Estimate of the biggest x at the solution. (default 1/N)
%    zsize      Estimate of the biggest z at the solution. (default 1/N)
%               xsize,zsize are used to scale (x,y,z).  Good estimates
%               should improve the performance of the barrier method
%
%  optParam     Structure with optimization parameters.
%  ========     The following fields are used:
% 
%    MaxIter    Maximum number of iterations (Prob.SOL.pdco.MaxIter)
%    MinorIter  Maximum number of iterations in LSQR (Prob.SOL.pdco.LSQRMaxIter)
%    eps_x      Accuracy for satisfying x1.*z1 = 0, x2.*z1 = 0,
%               where z = z1 - z2 and z1, z2 > 0.  (Prob.SOL.pdco.OptTol)
%    bTol       Accuracy for satisfying Ax + D2 r = b, A'y + z = gobj and x - x1
%               = bl, x + x2 = bu, where x1, x2 > 0. (Prob.SOL.pdco.FeaTol)
%    wait       = 0 solve the problem with default internal parameters;
%               = 1 pause, allows interactive resetting of parameters.
%               (Prob.SOL.pdco.wait)
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
%               The full [z;y] vector as returned from pdsco, including slacks
%               and extra linear constraints after rewriting constraints:
%               -inf < b_L < A*x < b_U < inf;  non-inf lower AND upper bounds.
%
%   ExitFlag    Tomlab Exit status from pdsco MEX.
%   Inform      pdsco information parameter: 0 = Solution found; 
%               1 = Too many iterations; 2 = Linesearch failed too often.
%
%   Iter        Number of iterations.
%   FuncEv      Number of function evaluations.
%   GradEv      Number of gradient evaluations.
%   HessEv      Number of Hessian evaluations.
%
%   Solver           Name of the solver (pdsco).
%   SolverAlgorithm  Description of the solver.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Jan 20, 2003.   Last modified Aug 17, 2006.

function Result=pdscoTL(Prob)

%#function pd_fgh

if nargin < 1, error('pdscoTL needs the Prob structure as input');return;end

Prob.solvType=3;

Prob = iniSolve(Prob,3,2,0);

Result=ResultDef(Prob);
Result.Solver = 'pdsco';

%LargeScale = DefPar(Prob,'LargeScale',1);

PriLev=DefPar(Prob,'PriLevOpt',0);

[bl,bu,n,nnLin,nnCon] = defblbu(Prob,inf,1);

if nnCon > 0
   fprintf('ERROR - PDSCO only solves linearly constrained problems');
   Result.ExitText = 'PDSCO does not allow nonlinear constraints';
   Result.ExitFlag = 10;
   Result=endSolve(Prob,Result);
   return
end

if isempty(Prob.HessPattern)
   % Only the diagonal of the Hessian is needed
   Prob.HessPattern = speye(n);
end

x_L = bl(1:n);
x_U = bu(1:n);
if any(x_L ~= 0) 
   fprintf('ERROR - PDSCO only solves x >= 0 problems\n');
   fprintf('All lower bounds x_L must be 0\n');
   Result.ExitText = 'PDSCO only solves x >= 0 problems';
   Result.ExitFlag = 10;
   Result=endSolve(Prob,Result);
   return
end
if any(~isinf(x_U)) 
   fprintf('ERROR - PDSCO only solves x >= 0 problems\n');
   fprintf('All upper bounds x_U must be inf\n');
   Result.ExitText = 'PDSCO only solves x >= 0 problems';
   Result.ExitFlag = 10;
   Result=endSolve(Prob,Result);
   return
end

nm = 0;
if nnLin > 0
   b_L = bl(n+1:n+nnLin);
   b_U = bu(n+1:n+nnLin);

   ixI = b_L ~= b_U;
   if any(ixI)
      ixL = ixI & ~isinf(b_L(1:nnLin));
      ixU = ixI & ~isinf(b_U(1:nnLin));
      ixB = ixL & ixU;
      jxL = find(ixL & ~ixB);
      jxU = find(ixU);
      if any(ixB)
         jxB = find(ixB);
         [iA,jA,vA]=find([Prob.A;Prob.A(jxB,:)]);
         iA=iA(:); jA=jA(:); vA=vA(:);
         nn = length(jxL) + length(jxU) + length(jxB);
         nm =  length(jxB);
         A = sparse([iA;jxU;jxL;nnLin+[1:nm]'],...
                    [jA;n+[1:nn]'],...
                    [vA;ones(length(jxU),1);-ones(length(jxL)+nm,1)], ...
                     nnLin+nm,n+nn);
         b = [b_U;b_L(jxB)];
         b(jxL)=b_L(jxL);
      else
         nn = length(jxL) + length(jxU);
         [iA,jA,vA]=find(Prob.A);
         iA=iA(:); jA=jA(:); vA=vA(:);
         A = sparse([iA;jxU;jxL],...
                    [jA;n+[1:nn]'],...
                    [vA;ones(length(jxU),1);-ones(length(jxL),1)], ...
                     nnLin,n+nn);
         b = b_U;
         b(jxL)=b_L(jxL);
      end
      x_L =  [x_L; zeros(nn,1)];
      x_U =  [x_U; inf*ones(nn,1)];
   else
      A = Prob.A;
      b = b_U;
      nn = 0;
   end
else
   %A = zeros(0,n); b_L = []; b_U = []; nn = 0; b = []; CRASH, must return
   % Instead add dummy constraint
   A   = sparse([zeros(1,n),1]); b_L = []; b_U = []; nn = 1; b = 1;
   nm  = 1;
   x_L =  [x_L; 0];
   x_U =  [x_U; inf];
end

N = n + nn;
M = nnLin + nm;

xsize = DefPar(Prob.SOL,'xsize',[]);

XS = 1;
% XS = 1/N;

% Starting point
if isempty(Prob.x_0)
   x_0 = ones(N,1)/N;   
   if isempty(xsize), xsize = XS; end
else
   x_0  = DefPar(Prob,'x_0',ones(n,1)/N);
   if nn > 0
      x_0  = [x_0;ones(nn,1)/N];
   end
   if isempty(xsize), xsize = max(abs(x_0)); end
end

% Safe-guard starting point
x_0   = max(x_L(1:N),min(x_0,x_U(1:N)));

if xsize == 0, xsize = XS; end

y0 = DefPar(Prob.SOL,'y0',[]);

if isempty(y0) 
   y0 = zeros(length(b),1);
end
if length(y0) < nnLin
   y0 = [y0;zeros(nnLin-length(y0),1)];
end
if length(y0) < length(b)
   y0 = [y0;-y0(jxB)];
end

z0 = DefPar(Prob.SOL,'z0',[]);
zsize = DefPar(Prob.SOL,'zsize',[]);
ZS = 1;
%ZS = 1/N;

if isempty(z0)
   z0 = ones(N,1)/N;   
   if isempty(zsize), zsize = ZS; end
end
if length(z0) < N
   z0 = [z0;ones(N-length(z0),1)];
end
%if isempty(zsize), zsize = max(max(z0),1/N); end
%if zsize == 0, zsize = 1/N; end
if isempty(zsize), zsize = max(max(z0),1); end
if zsize == 0, zsize = ZS; end

DefOpt = pdscoSet;

options = DefOpt;

UserSet = isfield(Prob.SOL,'pdco');
if UserSet
   UserOpt = Prob.SOL.pdco;
end

MaxIter = DefPar(Prob.optParam,'MaxIter',DefOpt.MaxIter);

if MaxIter ~= DefOpt.MaxIter
   options.MaxIter = MaxIter;
elseif UserSet
   options.MaxIter = DefPar(UserOpt,'MaxIter',DefOpt.MaxIter);
end

MinorIter = DefPar(Prob.optParam,'MinorIter',DefOpt.LSQRMaxIter);

if MinorIter ~= DefOpt.LSQRMaxIter
   options.LSQRMaxIter = MinorIter;
elseif UserSet
   options.LSQRMaxIter = DefPar(UserOpt,'LSQRMaxIter',DefOpt.LSQRMaxIter);
end

eps_x = DefPar(Prob.optParam,'eps_x',DefOpt.OptTol);

if eps_x ~= DefOpt.OptTol
   options.OptTol = eps_x;
elseif UserSet
   options.OptTol = DefPar(UserOpt,'OptTol',DefOpt.OptTol);
end

bTol = DefPar(Prob.optParam,'bTol',DefOpt.FeaTol);

if bTol ~= DefOpt.FeaTol
   options.FeaTol = bTol;
elseif UserSet
   options.FeaTol = DefPar(UserOpt,'FeaTol',DefOpt.FeaTol);
end

wait = DefPar(Prob.optParam,'wait',DefOpt.wait);

if wait ~= DefOpt.wait
   options.wait = wait;
elseif UserSet
   options.wait = DefPar(UserOpt,'wait',DefOpt.wait);
end

if UserSet
   options.gamma      = DefPar(UserOpt,'gamma',DefOpt.gamma);
   options.delta      = DefPar(UserOpt,'delta',DefOpt.delta);
   options.StepTol    = DefPar(UserOpt,'StepTol',DefOpt.StepTol);
  %options.StepSame   = DefPar(UserOpt,'StepSame',DefOpt.StepSame);
   options.x0min      = DefPar(UserOpt,'x0min',DefOpt.x0min);
   options.z0min      = DefPar(UserOpt,'z0min',DefOpt.z0min);
   options.mu0        = DefPar(UserOpt,'mu0',DefOpt.mu0);
   options.Method     = DefPar(UserOpt,'Method',DefOpt.Method);
   options.LSproblem  = DefPar(UserOpt,'LSproblem',DefOpt.LSproblem);
   options.LSQRatol1  = DefPar(UserOpt,'LSQRatol1',DefOpt.LSQRatol1);
   options.LSQRatol2  = DefPar(UserOpt,'LSQRatol2',DefOpt.LSQRatol2);
   options.LSQRconlim = DefPar(UserOpt,'LSQRconlim',DefOpt.LSQRconlim);
end

%options.mu0       = 1e-1;  % 1.0 starts near central path
%options.LSQRatol1 = 1e-6;  % Let LSQR solve loosely to start with

% Not used by pdsco:
%D1 = DefPar(Prob.SOL,'d1',1E-4);
%D2 = DefPar(Prob.SOL,'d2',1E-4);

Result.f_0  = nlp_f(x_0(1:n),Prob);

switch options.Method
   case 1
      Result.SolverAlgorithm = 'TOMLAB PDSCO, using Cholesky';
   case 2
      Result.SolverAlgorithm = 'TOMLAB PDSCO, using QR';
   case 3
      Result.SolverAlgorithm = 'TOMLAB PDSCO, using SOL Matlab LSQR';
  %case 5
  %   Result.SolverAlgorithm = 'TOMLAB PDSCO, Tlsqr with callbacks';
  %case 6
  %   Result.SolverAlgorithm = 'TOMLAB PDSCO, Iterative Tlsqr with callbacks';
   otherwise  % Case 4 default
      Result.SolverAlgorithm = 'TOMLAB PDSCO, using Tomlab Tlsqr Special Mex';
end

[x,y,z,inform,PDitns,CGitns,time] = ...
   pdsco( 'pd_fgh',A,b,M,N,options,x_0,y0,z0,xsize,zsize, Prob );

%plot(b-A*x)
%sum(max(0,b-A*x))

Result.Inform    = inform;
Result.Iter      = PDitns;
Result.MinorIter = CGitns;

Result.x_k  = x(1:n);
Result.f_k  = nlp_f(x(1:n),Prob);
Result.x_0  = x_0(1:n);
Result.y_k  = [z;y];

if ~isempty(Prob.FUNCS.g),Result.g_k  = nlp_g(x(1:n),Prob);end
if ~isempty(Prob.FUNCS.H),Result.H_k  = nlp_H(x(1:n),Prob);end

if nm > 0 & nnLin > 0
   v      = y(1:nnLin);
   v(jxB) = y(jxB) - y(nnLin+1:nnLin+nm); % Possible sign error?
   Result.v_k  = [z(1:n);v];
else
   Result.v_k  = [z(1:n);y];
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

% Must check x and Ax with looser tolerance
%eps_x = DefPar(options,'OptTol',1.73E-6);
%bTol = DefPar(options,'LSQRatol',1E-8);

Result = StateDef(Result, x(1:n), Result.Ax, [], eps_x, bTol, [], bl, bu,1);

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
% 030818 ango Fixed bug when Prob.A is one row
% 031204 ango Changed pd_fgH reference to pd_fgh
% 040110 hkh  Revised for v4.2
% 040809 med  Pragmas added for MATLAB Compiler
% 041201 hkh  Call endSolve before error return
% 041202 hkh  Revise calls to defblbu and StateDef, use output from defblbu
% 041222 med  Safeguard added for x_0
% 050109 med  Safeguard corrected
% 050726 med  x_0 safeguard corrected again
% 060814 med  FUNCS used for callbacks instead
% 060817 hkh  PDSCO needs linear constraints, add dummy constraint if none
