% TOMLAB filterSQP nonlinear constrained optimization solver
%
% function Result = filterSQPTL(Prob)
%
% filterSQP solves constrained nonlinear constrained problems of
% the type:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds
%    b_L <= A*x  <= b_U, linear constraints
%    c_L <= c(x) <= c_L, nonlinear constraints
%
%
% ---------------------------------------------------------------------
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%  x_L, x_U     Bounds on variables.
%  b_L, b_U     Bounds on linear constraints.
%  c_L, c_U     Bounds on nonlinear constraints.
%               For equality constraints (or fixed variables), set
%               e.g. b_L(k) == b_U(k)
%
%  LargeScale   If 1 use sparse version of solver. The default is 0, the
%               dense version.
%
%  PriLevOpt    Print level in the filterSQP solver.
%
%  WarmStart    Indicates that the solver should be
%               warmstarted. See Prob.DUNDEE for necessary
%               arguments when doing warmstarts.
%
%
%  optParam     Structure with optimization parameters. Fields
%               used:
%
%    MaxIter    Maximum number of iterations.
%
%
%  DUNDEE       Structure with special fields for filterSQP optimization
%               parameters. The following fields are used:
%
%    QPmin      Lower bound for the QP subproblems.
%               If not set, use Prob.f_Low. Default: -1E300
%
%    rho        Initial trust region radius.
%               Default: 10.0 (REAL)
%
%    kmax       Maximum size of the null-space, less than or equal to no. of variables
%               Default: n    (INTEGER)
%
%    maxf       Maximum size of the filter
%               Default: 100  (INTEGER)
%
%    mlp        Maximum level parameter for resolving degeneracy in BQPD QP subsolver.
%               Default: 100  (INTEGER)
%
%    Name       Problem name, at most 10 characters. The output files
%               are named <pname>.sum and <pname>.out
%               Default name filterSQP, i.e. files filterSQP.sum, filterSQP.out
%
%
%    optPar     Vector of max length 20 with optimization parameters:
%               If any element is -999, default value is assigned.
%
%
%     optPar(1):  iprint  0     Print level in filterSQP
%                               Summary on file filterSQP.sum (if default)
%                               More printout on file filterSQP.out (if default)
%     optPar(2):  tol     1E-10 Relative tolerance for BQPD subsolver
%     optPar(3):  emin    1.0   1=Use cscale in BQPD, 0=no scaling
%     optPar(4):  sgnf    5E-4  Max rel error in two numbers equal in
%                              exact arithmetic (BQPD)
%
%     optPar(5):  nrep    2     Max number of refinement steps (BQPD)
%     optPar(6):  npiv    3     No repeat if no more than npiv steps were taken
%     optPar(7):  nres    2     Max number of restarts if unsuccessful
%     optPar(8):  nfreq   500   The max interval between refactorizations
%
%     optPar(9):  NLP_eps 1E-6  NLP subproblem tolerance
%     optPar(10): ubd     1E-2  Upper bound on constraint violation used
%                               in the filter
%     optPar(11): tt      0.125 Parameter related to ubd. The actual
%                               upper bound is defined by the maximum of
%                               ubd and tt multiplied by the initial
%                               constraint violation
%
%     optPar(19): infty   1E20  A large value representing infinity
%
%     optPar(20): Nonlin  0     If 1, skip linear feasibility tests
%                               filterSQP treating all constraints as nonlinear
%
%
%    lws        If doing warmstarts, this field is set to the
%               Result.DUNDEE.lws field from the previous run.
%
%    istat      Similarly, for warmstarts, set istat to
%               Result.DUNDEE.istat from the previous run. Only the
%               first element is used.
%
%    lam        Vector of initial multipliers. Necessary for
%               warmstarts, but can always be given if desired.
%               Must be n+m elements in order to be used.
%
%    morereal   Increase of REAL workspace. A problem dependent
%               default value is used if <0 or empty.
%
%    moreint    Increase of INTEGER workspace. A problem dependent
%               default value is used if <0 or empty.
%
%      Scaling parameters: It is possible to supply scale factors for
%      the variables and/or the constraints. Normally, the DUNDEE
%      solvers does not differentiate between linear and nonlinear
%      constraints with regard to scaling, but the Tomlab interface
%      handles this automatically. Thus is it possible to give scale
%      factors e.g. for the nonlinear constraints only. All scaling
%      values must be greater than zero.
%
%      The three parameters in the Prob.DUNDEE substructure that
%      control scaling are:
%
%    xScale     Vector of scale factors for variables. If less than
%               n values given, 1's are used for the missing elements.
%
%    bScale     Vector of scale factors for the linear constraints. If
%               length(bScale) is less than the number of linear
%               constraints ( size(Prob.A,1) ), 1's are used for the
%               missing elements.
%
%    cScale     Vector of scale factors for the nonlinear
%               constraints. If length(cScale) is less than the
%               number of nonlinear constraints, 1's are used for
%               the missing elements.
%
% ----------------------------------------------------------------------------------
%
% OUTPUT:
% Result        Structure with optimization results
%
%   f_k         Function value at optimum.
%   g_k         Gradient of the function.
%
%   x_k         Solution vector.
%   x_0         Initial solution vector.
%
%   c_k         Nonlinear constraint residuals.
%   cJac        Nonlinear constraint gradients.
%
%
%   xState      State of variables. Free == 0; On lower == 1; On upper == 2;
%               Fixed == 3;
%
%   bState      State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   cState      State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   v_k         Lagrangian multipliers (for simple bounds +
%               dual solution vector (linear and nonlinear constraints).
%
%   ExitFlag    Exit status from filterSQP MEX.
%   Inform      filterSQP information parameter:
%                0 - Solution found
%                1 - Unbounded problem. Feasible x with f(x)<=fmin found
%                2 - Linear constraints are inconsistent
%                3 - (Locally) nonlinear infeasible, optimal solution to feasibility problem found
%                4 - Terminated at point with h(x)<=eps but QP infeasible
%                5 - Terminated with rho < eps
%                6 - Too many iterations
%                7 - Crash in user routine could not be resolved
%                8 - Unexpected ifail from QP solver
%                    This is often due to too little memory being
%                    allocated and is remedied by setting
%                    appropriate values in the Prob.DUNDEE.morereal
%                    and Prob.DUNDEE.moreint parameters.
%
%                9 - Not enough REAL workspace
%               10 - Not enough INTEGER workspace
%
%   Iter        Number of iterations.
%   FuncEv      Number of function evaluations.
%   GradEv      Number of gradient evaluations.
%   ConstrEv    Number of constraint evaluations.
%
%   Solver           Name of the solver (filterSQP).
%   SolverAlgorithm  Description of the solver.
%
%
% Result.DUNDEE   Substructure with filterSQP specific data
%
%   istat       Solution statistics, integer values. First element is
%               required as input if doing a warmstart.
%
%   lws         Workspace vector, should be treated as integer
%               valued. Required if doing warmstarts.
%
%   lam         Vector of multipliers, required if doing warmstarts.
%
%   istat(1)    Dimension of nullspace at solution
%   istat(2)    Number of iterations
%   istat(3)    Number of feasibility iterations
%   istat(4)    Number of objective evaluations
%   istat(5)    Number of constraint evaluations
%   istat(6)    Number of gradient evaluations
%   istat(7)    Number of Hessian evaluations
%   istat(8)    Number of QPs with mode<=2
%   istat(9)    Number of QPs with mode>=4
%   istat(10)   Total number of QP pivots
%   istat(11)   Number of SOC steps
%   istat(12)   Maximum size of filter
%   istat(13)   Maximum size of Phase 1 filter
%   istat(14)   Number of QP crashes
%
%   rstat      Solution statistics, floating point values.
%
%   rstat(1)    l_2 norm of KT residual
%   rstat(2)
%   rstat(3)    Largest modulus multiplier
%   rstat(4)    l_inf norm of final step
%   rstat(5)    Final constraint violation h(x)
%   rstat(6)
%   rstat(7)

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Dec 13, 2002.   Last modified May 26, 2007.

function Result=filterSQPTL(Prob)

%#function nlp_dcD nlp_d2L nlp_dcDS

if nargin < 1, error('filterSQPTL needs the Prob structure as input');return;end

Prob.solvType = 3; % NLP (CON) solver

if abs(Prob.NumDiff)>10
   Prob = iniSolve(Prob,3,2,1);
else
   Prob = iniSolve(Prob,3,1,1);
end

Result=ResultDef(Prob);
Result.Solver = 'filterSQP';

LargeScale = DefPar(Prob,'LargeScale',0);

switch(LargeScale)
case 1,
   Result.SolverAlgorithm = 'Sparse Filter SQP Trust Region';
otherwise,
   Result.SolverAlgorithm = 'Dense Filter SQP Trust Region';
end

PriLev=DefPar(Prob,'PriLevOpt',0);

optPar = DefPar(Prob.DUNDEE,'optPar',-999*ones(20,1));
if(length(optPar)<20), optPar(end+1:20)=-999; end

% Set infty (BIG) = Real value for infinity, to default 1E20
if(optPar(19)<0), optPar(19) = 1E20; end

% Safe guard BIG, always >= 1E20
optPar(19) = max(1E20,optPar(19));
BIG        = optPar(19);

% Get single bounds vectors
% defblbu Order 0 gives [ bounds nonlinear linear]
[bl,bu,n,nnLin,nnCon] = defblbu(Prob,BIG,0);

m = nnCon+nnLin;

% Starting point
x_0        = DefPar(Prob,'x_0',zeros(n,1));

% Safe guard x_0 
x_0 = max(bl(1:n),min(bu(1:n),x_0));

Result.f_0 = nlp_f(x_0, Prob);
Result.x_0 = x_0;

% Setup augmented a matrix

if LargeScale
   if isempty(Prob.ConsPattern)
      a = sparse( [ones(1,n) ; ones(nnCon,n) ; Prob.A]' );
   else
      a = sparse( [ones(1,n) ; Prob.ConsPattern ; Prob.A]' );
   end
else
   a = full( [ones(1,n) ; ones(nnCon,n) ; Prob.A]' );
end

% if nnCon > 0
%    % Determine the sparse problem structure
%    if ~isempty(Prob.ConsPattern)
%       [ix,iy]=find(Prob.ConsPattern');
%       
%       % Send linear index from multiple subscripts for nonzero pattern
%       Prob.ConsIdx = sub2ind(size(Prob.ConsPattern'),ix,iy);
%    end
% end

% Warmstart. Send istat(1), lws and lam back into the MEX.
WarmStart = DefPar(Prob,'WarmStart',double(0));

if WarmStart==1
  istat = DefPar(Prob.DUNDEE,'istat',[]);
  if isempty(istat)
    error('filterSQPTL: Warmstart requested but Prob.DUNDEE.istat is empty')
  else
    % Just first element is required
    istat = istat(1);
  end
  
  lws = DefPar(Prob.DUNDEE,'lws',[]);
  if isempty(lws)
    error('filterSQPTL: Warmstart requested but Prob.DUNDEE.lws is empty');
  end

  lam = DefPar(Prob.DUNDEE,'lam',[]);
  if isempty(lam) 
    error(['filterSQPTL: Warmstart requested but Prob.DUNDEE.lam is' ...
	   ' empty']);
  end

else
  lam   = DefPar(Prob.DUNDEE,'lam',zeros(n+m,1));
  lws   = [];
  istat = [];
end

% Scaling parameters

% First assume no scaling at all. DUNDEE solvers only provide
% "constraint scaling", not differentiating between linear and
% nonlinear scaling factors, so if at least one of the linear or
% nonlinear scale factors are given, we need to scale BOTH types of
% constraints.

smode = zeros(3,1);

xScale = DefPar(Prob.DUNDEE,'xScale',[]);
bScale = DefPar(Prob.DUNDEE,'bScale',[]);
cScale = DefPar(Prob.DUNDEE,'cScale',[]);

xScale = double(xScale(:)); 
bScale = double(bScale(:)); 
cScale = double(cScale(:)); 

% Check each of the scaling vectors for values < 0
% These can be fixed by taking absolute values
i=find(xScale<0);
if any(i)
   if(PriLev>0)
      warning(['Negative elements found in DUNDEE.xScale. Absolute' ...
            ' values are substituted.']);
   end
   xScale(i) = abs(xScale(i));
end
i=find(bScale<0);
if any(i)
   if(PriLev>0)
      warning(['Negative elements found in DUNDEE.bScale. Absolute' ...
            ' values are substituted.']);
   end
   bScale(i) = abs(bScale(i));
end
i=find(cScale<0);
if any(i)
   if(PriLev>0)
      warning(['Negative elements found in DUNDEE.cScale. Absolute' ...
            ' values are substituted.']);
   end
   cScale(i) = abs(cScale(i));
end

% Zeros are worse - error if any zeros
if(any(xScale==0) | any(bScale==0) | any(cScale==0))
   error('filterSQPTL: one or more scale factors are zero');
end

% Check each scale factor vector, padding with zeros if needed.
if ~isempty(xScale)
  smode(1) = 1;
  
  j = length(xScale);
  if j<n 
    % Pad with ones if too few factors given
    xScale = [ xScale ; ones(n-j,1) ];
  else
    % Use no more than n factors
    xScale = xScale(1:n);
  end

  % If all scale factors are one, there is no need to scale
  if( all(xScale==1.0) )
    smode(1)=0;
  else
    smode(1)=1;
  end
  
else
  % Nothing given - use ones, but disable variable scaling
  smode(1) = 0;
  xScale = ones(n,1);
end

% Nonlinear constraints scaling
if ~isempty(cScale)
  j = length(cScale);
  if j<nnCon
    % Pad with ones if less than nnCon factors given
    cScale = [ cScale ; ones(nnCon-j,1) ];
  else
    % No more than nnCon elements used
    cScale = cScale(1:nnCon);
  end
  
  % If all scale factors are one, there is no need to scale
  if( all(cScale==1.0) )
    smode(2)=0;
  else
    smode(2)=1;
  end
  
else
  % Nothing given - use ones, but indicate that (so far) constraint
  % scaling is <OFF>.
  smode(2) = 0;
  cScale = ones(nnCon,1);
end

  
% Linear constraints scaling
if ~isempty(bScale)
  j = length(bScale);
  if j<nnLin
    % Pad with ones if less than nnLin factors given
    bScale = [bScale ; ones(nnLin-j,1)];
  else
    % At most nnLin elements are used
    bScale = bScale(1:nnLin);
  end
  
  % If all factors are one, there is no need to scale
  if( all(bScale==1.0) )
    smode(3)=0;
  else
    smode(3)=1;
  end
else
  % Nothing given - use ones. Set smode(3)=0 to indicate that we
  % may not need to do any scaling of constraints.
  smode(3) = 0;
  bScale   = ones(nnLin,1);
end

% Now put together a scalar parameter telling the scaling mode:
% 0 - no scaling
% 1 - variable scaling
% 2 - constraint scaling only
% 3 - variable AND constraint scaling

scmode = 0;
if smode(1)>0, scmode=1; end
if (smode(2)>0 | smode(3)>0), scmode=scmode+2; end

% All the scale factors together
Scale = double([xScale;cScale;bScale]);

fLow = DefPar(Prob.DUNDEE,'QPmin',Prob.f_Low);
kmax = min(n,DefPar(Prob.DUNDEE,'kmax',n ));
mlp  = max(3,min(m+1,DefPar(Prob.DUNDEE,'mlp',100)));
maxf = DefPar(Prob.DUNDEE,'maxf',100);
rho  = DefPar(Prob.DUNDEE,'rho',10);

MaxIter = DefPar(Prob.optParam,'MaxIter',2000);

% Base name for print and summary files
pname = DefPar(Prob.DUNDEE,'Name','filterSQP');
if length(pname)>10, pname=pname(1:10); end


% moremem parameter
moremem(1) = DefPar(Prob.DUNDEE,'morereal',-1);
moremem(2) = DefPar(Prob.DUNDEE,'moreint',-1);

if LargeScale == 1
   [ifail,x_k,f_k,c_k,v_k,lws,istat,rstat] = ...
      filSQPs(a,bl,bu,nnCon,x_0,Scale,scmode,fLow,MaxIter,...
      rho,mlp,kmax,maxf,WarmStart,lws,lam,istat,PriLev,pname,...
      optPar,Prob,moremem);
else
   [ifail,x_k,f_k,c_k,v_k,lws,istat,rstat] = ...
      filSQPd(a,bl,bu,nnCon,x_0,Scale,scmode,fLow,MaxIter,...
      rho,mlp,kmax,maxf,WarmStart,lws,lam,istat,PriLev,pname,...
      optPar,Prob,moremem);
end


Result.Inform = ifail;

Result.x_k  = x_k;
Result.f_k  = f_k;
Result.x_0  = x_0;

if ~isempty(Prob.FUNCS.g),Result.g_k  = nlp_g(x_k,Prob);end
if ~isempty(Prob.FUNCS.H),Result.H_k  = nlp_H(x_k,Prob);end

Result.v_k  = v_k;

Result.DUNDEE.istat = istat;
Result.DUNDEE.rstat = rstat;
Result.DUNDEE.ifail = ifail;
Result.DUNDEE.lws   = lws;
Result.DUNDEE.lam   = v_k;

Result.Iter     = istat(2);
% Result.FuncEv   = istat(4);
% Result.ConstrEv = istat(5);
% Result.GradEv   = istat(6);
% Result.HessEv   = istat(7);

optParam = Prob.optParam;

if(nnLin>0)
  %Ax = Prob.A*x_k;
  Result.Ax = c_k(nnCon+1:m);
else
  Result.Ax = [];
end

if(nnCon>0)
  Result.c_k = c_k(1:nnCon);
  if ~isempty(Prob.FUNCS.dc)
     Result.cJac = nlp_dc(x_k,Prob);
  end
end

Result = StateDef(Result, x_k, Result.Ax, Result.c_k, optParam.xTol, ...
                  optParam.bTol, optParam.cTol, bl, bu, 0);

switch Result.Inform
   case 6
     ExitFlag=1;  % Too many iterations
   case 1
     ExitFlag=2;  % Unbounded
   case {2,3,4}
     ExitFlag=4;  % Infeasible
   case {5,8}
     ExitFlag=3;  % Rank problem
   case {9,10,7}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end
Result.ExitFlag = ExitFlag;

switch(Result.Inform)
case 0,
   Result.ExitText = 'Solution found';
case 1,
   Result.ExitText = 'Unbounded problem. Feasible x with f(x)<=fmin found';
case 2,
   Result.ExitText = 'Linear constraints are inconsistent';
case 3,
   Result.ExitText = '(Locally) nonlinear infeasible, optimal solution to feasibility problem found';
case 4,
   Result.ExitText = 'Terminated at point with h(x)<=eps but QP infeasible';
case 5,
   Result.ExitText = 'Terminated with rho < eps';
case 6,
   Result.ExitText = 'Too many iterations';
case 7,
   Result.ExitText = 'Crash in user routine could not be resolved';
case 8,
   Result.ExitText = 'Unexpected ifail from QP solver';
case 9,
   Result.ExitText = 'Not enough REAL workspace';
case 10,
   Result.ExitText = 'Not enough INTEGER workspace';
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 021223 hkh  Avoid calls to Prob.FUNCS.f etc. directly, could crash
% 021223 hkh  Use Prob.ConsPattern, no local ConsPattern
% 021223 hkh  Default MaxIter 1000, optPar(1:end) = -999 as default
% 021223 hkh  Compute pattern based on transpose of ConsPattern
% 021223 hkh  Change computation of Result.ExitFlag, Result.c_k, Result.Ax
% 030103 ango Change scaling parameters
% 030113 ango Add parameter DUNDEE.QPmin
% 030114 ango Add optPar(1) - iprint
% 030128 ango Revised comments
% 030128 hkh  Switching of bl and bu not correct, use nnLin, not nnCon
% 030206 hkh  Add optPar(20)=Nonlin. If 1, no special treatment of linear cons
% 030206 hkh  Change optPar comments, change name of file to filterSQPTL
% 030206 hkh  Set Result.x_0
% 030220 hkh  Restrict mlp and kmax,  min(m, mlp), min(n,kmax)
% 031114 hkh  Change infty,BIG to 1E12 instead of 1E20, set default optPar(19)
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040713 hkh  Change default MaxIter to 2000, corresponding to optParamDef
% 041202 hkh  Revise calls to defblbu and StateDef, use different Order=0
% 041202 hkh  filterSQP fails if BIG < 1D20, correct this bug, safe guard BIG
% 050829 med  Corrected FuncEv, GradEv, HessEv and ConstrEv (set by TOMLAB)
% 060714 med  Spelling updated
% 060814 med  FUNCS used for callbacks instead
% 060818 hkh  Use Prob.f_Low instead of -1E300 defining fLow, if QPmin not set 
% 060818 hkh  Result.v_k should be full v_k, not Result.v_k  = v_k(n+1:end);
% 060818 med  isnan checks removed for x_0
% 060829 med  LargeScale added in help
% 070526 med  Output removed