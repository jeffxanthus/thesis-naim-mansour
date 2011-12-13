% TOMLAB minlpBB MINLP Solver
%
% function Result = minlpBBTL(Prob)
%
% minlpBB solves constrained nonlinear mixed-integer problems of
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
%    x_i are restricted to integer values for i in I,
%        where I is a subset of {1,2,...,n}, x in R^n
%
%    In addition, Special Ordered Sets of type 1 (SOS1) can be
%    defined.
%
% The algorithm uses a branch-and-bound scheme with a depth-first
% search strategy. The NLP relaxations are solved using the solver
% filterSQP by R.Fletcher and S.Leyffer.
%
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%    A          Linear constraints coefficient matrix
%    x_L, x_U   Bounds on variables.
%    b_L, b_U   Bounds on linear constraints.
%    c_L, c_U   Bounds on nonlinear constraints.
%               For equality constraints (or fixed variables), set
%               e.g. b_L(k) == b_U(k).
%
%    LargeScale If 1 use sparse version of solver. The default is 0, the
%               dense version.
%
%    PriLevOpt  Print level in MEX interface and minlpBB solver.
%
% Prob.optParam Structure with optimization parameters
%
%    MaxIter    Limit of iterations
%
% Prob.MIP      Structure with fields defining the integer properties of
%               the problem. The following fields are used:
%
%   IntVars:
%               If empty, all variables are assumed non-integer
%               If islogical(IntVars) (=all elements are 0/1), then
%               1 = integer variable, 0 = continuous variable.
%               If any element >1, IntVars is the indices for integer variables
%
%   VarWeight   Defines the priorities of the integer variables.
%               Can be any values, but minlpBB uses integer
%               priorities internally, with higher values implying
%               higher priorities.
%
%
%   sos1        Structure defining the Special Ordered Sets of Type 1 (SOS1).
%               If there are k sets of type sos1, then
%               sos1(1).var is a vector of indices for variables in sos1, set 1.
%               sos1(1).row is the row number for the reference row identifying
%               the ordering information for the sos1 set, i.e.
%               A(sos1(1).row,sos1(1).var) identifies this information
%               sos1(1).prio sets the priority for sos1 test 1.
%
%               sos1(2).var is a vector of indices for variables in sos1, set 2.
%               sos1(2).row is the row number for the reference row of sos1 set 2.
%               sos1(2).prio is the priority for sos1 set 2.
%                  ...
%               sos1(k).var is a vector of indices for variables in sos1, set k.
%               sos1(k).row is the row number for the reference row of sos1 set k.
%               sos1(k).prio is the priority for sos1 set k.
%
%
% Prob.DUNDEE   Structure with special fields for minlpBB optimization parameters.
%               The following fields are used:
%
%    stackmax   Maximum size of the LIFO stack storing info about B&B tree.
%               Default 10000
%
%    QPmin      Lower bound for the QP subproblems.
%               If not set, use Prob.f_Low. Default: -1E300
%
%    rho        Initial trust region radius.
%               Default: 10.0 (REAL)
%
%    kmax       Maximum size of the null-space, less than or equal to no. of variables.
%               Default: n    (INTEGER)
%
%    maxf       Maximum size of the filter.
%               Default: 100  (INTEGER)
%
%    mlp        Maximum level parameter for resolving degeneracy in BQPD QP subsolver.
%
%    lam        Multipliers (n+m) on entry (NOTE: Experimental parameter).
%
%    Name       Problem name, at most 10 characters. The output files
%               are named <pname>.sum and <pname>.out.
%               Default name minlpBB, i.e. files minlpBB.sum, minlpBB.out
%
%
%   optPar      Vector of max length 20 with optimization parameters.
%               If any element is -999, default value is assigned.
%               The elements used by minlpBB are:
%
%     optPar(1):  iprint  0     Print level in minlpBB
%                               Summary on file minlpBB.sum
%                               More printout on file minlpBB.out
%     optPar(2):  tol     1E-10 Relative tolerance for BQPD subsolver
%     optPar(3):  emin    1.0   1=Use cscale in BQPD, 0=no scaling
%     optPar(4):  sgnf    5E-4  Max rel error in two numbers equal in
%                               exact arithmetic (BQPD)
%     optPar(5):  nrep    2     Max number of refinement steps (BQPD)
%     optPar(6):  npiv    3     No repeat if no more than npiv steps were taken
%     optPar(7):  nres    2     Max number of restarts if unsuccessful
%     optPar(8):  nfreq   500   The max interval between refactorizations
%
%     optPar(9):  ubd     1E-2  Upper bound on constraint violation used
%                               in the filter
%     optPar(10): tt      0.125 Parameter related to ubd. The actual
%                               upper bound is defined by the maximum of
%                               ubd and tt multiplied by the initial
%                               constraint violation
%     optPar(11): NLP_eps  1E-6 NLP subproblem tolerance
%     optPar(12): epsilon  1D-6 Tolerance for x-value tests,
%                               and constraint feasibility
%     optPar(13): MIopttol 1E-4 Tolerance for function value tests
%
%     optPar(17): branchtype    Branch strategy. Currently only 1 strategy
%     optPar(19): infty    1E20 A large value representing infinity
%                               Should always be >= 1E20
%
%     optPar(20): Nonlin 0      If 1, skip linear feasibility tests
%                               minlpBB treating all constraints as nonlinear
%
%
%    morereal   Number of extra REAL workspace locations. Set to <0
%               for problem dependent default strategy.
%
%    moreint    Number of extra INTEGER workspace locations. Set to <0
%               for problem dependent default strategy.
%
%      Scaling parameters. It is possible to supply scale factors for
%      the variables and/or the constraints. Normally, the DUNDEE
%      solvers does not differentiate between linear and nonlinear
%      constraints with regard to scaling, but the TOMLAB interface
%      handles this automatically. Thus is it possible to give scale
%      factors e.g. for the nonlinear constraints only. The three
%      parameters in the Prob.DUNDEE substructure that control scaling are:
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
% OUTPUT:
% Result   Structure with optimization results
%
%   f_k      Function value at optimum.
%   g_k      Gradient of the function.
%
%   x_k      Solution vector.
%   x_0      Initial solution vector.
%
%   c_k      Nonlinear constraint residuals.
%   cJac     Nonlinear constraint gradients.
%
%   xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%            Fixed == 3;
%
%   bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%            Equality == 3;
%
%   cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%            Equality == 3;
%
%   v_k        Lagrangian multipliers (for bounds + dual solution vector).
%
%   ExitFlag   Exit status.
%   Inform     minlpBB information parameter
%               0 - Optimal solution found
%               1 - Root problem infeasible
%               2 - Integer infeasible
%               3 - Stack overflow - some integer solution. obtained
%               4 - Stack overflow - no integer solution obtained
%               5 - SQP termination with rho < eps
%               6 - SQP termination with iter > max_iter
%               7 - Crash in user supplied routines
%               8 - Unexpected ifail from QP solvers
%                   This is often due to too little memory being
%                   allocated and is remedied by setting
%                   appropriate values in the Prob.DUNDEE.morereal
%                   and Prob.DUNDEE.moreint parameters.
%
%               9 - Not enough REAL workspace or parameter error
%              10 - Not enough INTEGER workspace or parameter error
%
%
%   rc               Reduced costs. If ninf=0, last m == -v_k.
%   Iter             Number of iterations.
%   FuncEv           Number of function evaluations.
%   GradEv           Number of gradient evaluations.
%   ConstrEv         Number of constraint evaluations.
%   QP.B             Basis vector in TOMLAB QP standard.
%   Solver           Name of the solver (minlpBB).
%   SolverAlgorithm  Description of the solver (sparse or dense, mainly).

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Dec 13, 2002.   Last modified May 26, 2007.

function Result=minlpBBTL(Prob)

%#function nlp_dcDS nlp_d2L nlp_dcD

if nargin < 1, error('minlpBBTL needs the Prob structure as input');end

Prob.solvType = 12; % MINLP solver

if abs(Prob.NumDiff)>10
   Prob = iniSolve(Prob,3,2,1);
else
   Prob = iniSolve(Prob,3,1,1);
end

Result=ResultDef(Prob);
Result.Solver = 'minlpBB';

LargeScale = DefPar(Prob,'LargeScale',0);

switch(LargeScale)
case 1,
   Result.SolverAlgorithm = 'Sparse Branch and Bound MINLP';
otherwise,
   Result.SolverAlgorithm = 'Dense Branch and Bound MINLP';
end

PriLev = DefPar(Prob,'PriLevOpt',0);

optPar = DefPar(Prob.DUNDEE,'optPar',-999*double(ones(20,1)));
if(length(optPar)<20), optPar(end+1:20)=-999; end

% Set infty (BIG) = Real value for infinity, to default 1E20
if(optPar(19)<0), optPar(19) = 1E20; end

% Safe guard BIG, always >= 1E20
optPar(19) = max(1E20,optPar(19));
BIG        = optPar(19);

[bl,bu,n,nnLin,nnCon] = defblbu(Prob,BIG,0);

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
      error('minlpBBTL: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

VW = DefPar(Prob.MIP,'VarWeight',[]);

if ~isempty(VW)
   [i1,i2]=sort(VW);
   Priority(i2(i2)) = -i2;
   Priority = Priority - min(Priority);
else
   Priority = zeros(size(IntVars));
end

% SOS1 sets
sos1 = DefPar(Prob.MIP,'sos1',[]);

if isempty(sos1)
   setbeg    = [];
   setcolidx = [];
   setref    = [];
   setrefrow = [];
   setprio   = [];
else
   ns1       = length(sos1);
   nset      = ns1;
   setbeg    = 1;
   setcolidx = [];
   setref    = [];
   setrefrow = [];
   setprio   = [];
   
   % index for removing rows in A and b_L, b_U
   Aidx=1:nnLin;
   
   for i=1:ns1
      if ~isfield(sos1(i),'var')
         fprintf('sos1 set %d. ',i);
         error('minlpBBTL: sos1 field var is missing');
      end
      if ~isfield(sos1(i),'row')
         fprintf('sos1 set %d. ',i);
         error('minlpBBTL: sos1 field row is missing');
      end
      
      ix = sos1(i).var;

      if ~(all(ix >= 1 & ix <= n))
         fprintf('sos1 set %d. ',i);
         fprintf('\n');
         error('minlpBBTL: Illegal sos1 input variable vector');
      end
      
      row = sos1(i).row;
      if ~(row >= 0 & row <= nnLin)
         fprintf('sos1 set %d. ',i);
         fprintf('Illegal row number  %d.',row);
         fprintf('\n');
         error('minlpBBTL: Illegal sos1 row data');
      end
      
      prio = DefPar(sos1(i),'prio',0);
      
      k         = length(ix);
      setbeg    = [setbeg; setbeg(length(setbeg))+k];
      setcolidx = [setcolidx; ix(:)];
      setprio   = [setprio; prio];
      
      if row==0
         refrow = full(Prob.QP.c(ix));
         setrefrow = [setrefrow; refrow(:)];
      else
         refrow = full(Prob.A(row,ix));
         setrefrow = [setrefrow; refrow(:)];
         
         % Remove this row from A
         Aidx(row) = 0;
      end
   end
   
   % Remove the appropriate rows from A and b_L, b_U OR NOT? Most
   % problems fail/do not solve correctly if rows are taken out
   %
   %     Aidx = Aidx( find(Aidx) );
   %     Prob.A   = Prob.A( Aidx,: );
   %     Prob.b_L = Prob.b_L( Aidx );
   %     Prob.b_U = Prob.b_U( Aidx );
   %     ix       = [1:n + nnCon;n+nnCon+Aidx];
   %     bl       = bl(ix);
   %     bu       = bu(ix);
   %     nnLin    = length(Aidx);
end

m = nnLin+nnCon;

% Starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));

% Safe guard x_0
x_0 = max(bl(1:n),min(bu(1:n),x_0)); 

Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

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
   error('minlpBBTL: one or more scale factors are zero');
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
if smode(1)>0, scmode = 1; end
if (smode(2)>0 | smode(3)>0), scmode = scmode+2; end

% All the scale factors together
Scale = double([xScale;cScale;bScale]);

% Previous scaling code was just these two lines:
% Scale  = DefPar(Prob.scale',double(ones(n+m,1)));
% scmode = DefPar(Prob.DUNDEE,'scmode',0);

fLow = DefPar(Prob.DUNDEE,'QPmin',Prob.f_Low);
kmax = min(n,DefPar(Prob.DUNDEE,'kmax',n ));

% Maximum number of levels of recursion, typically 20, maximum = m.
mlp  = max(3,min(m+1,DefPar(Prob.DUNDEE,'mlp',100)));

maxf = DefPar(Prob.DUNDEE,'maxf',100);

rho  = DefPar(Prob.DUNDEE,'rho',10);
smax = DefPar(Prob.DUNDEE,'stackmax',10000);

MaxIter   = DefPar(Prob.optParam,'MaxIter',5000);

% Memory increase parameters
moremem(1) = DefPar(Prob.DUNDEE,'morereal',-1);
moremem(2) = DefPar(Prob.DUNDEE,'moreint',-1);

% Multipliers
lam = DefPar(Prob.DUNDEE,'lam',zeros(n+m,1));

% Setup augmented a matrix
if LargeScale
   if isempty(Prob.ConsPattern)
      a = [spones(ones(1,n)) ; spones(ones(nnCon,n)) ; sparse(Prob.A)]';
   else
      a = [spones(ones(1,n)) ; sparse(Prob.ConsPattern) ; sparse(Prob.A)]';
   end
else
   a = [ones(1,n) ; ones(nnCon,n) ; full(Prob.A)]';
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

pname = DefPar(Prob.DUNDEE,'Name','minlpBB');
if length(pname) > 10, pname=pname(1:10); end

% Call the MEX file

if LargeScale == 1
  [ifail,x_k,f_k,c_k,v_k,iter_nlp] = minlpBBs(...
      a,bl,bu,nnCon,x_0,Scale,scmode,MaxIter,rho,...
      kmax,mlp,maxf,PriLev,pname,optPar,Prob,IntVars,Priority,...
      smax,setbeg,setcolidx,setrefrow,setprio,moremem,fLow,lam);
  
else
  [ifail,x_k,f_k,c_k,v_k,iter_nlp] = minlpBBd(...
      a,bl,bu,nnCon,x_0,Scale,scmode,MaxIter,rho,...
      kmax,mlp,maxf,PriLev,pname,optPar,Prob,IntVars,Priority,...
      smax,setbeg,setcolidx,setrefrow,setprio,moremem,fLow,lam);
end

Result.Inform = ifail;

Result.x_k  = x_k;
Result.f_k  = f_k;

if ~isempty(Prob.FUNCS.g), Result.g_k = nlp_g(x_k,Prob); else Result.g_k =[]; end
if ~isempty(Prob.FUNCS.H), Result.H_k = nlp_H(x_k,Prob); else Result.H_k =[]; end

if nnCon > 0
   Result.c_k = nlp_c(x_k,Prob); 
   if ~isempty(Prob.FUNCS.dc)
      Result.cJac = nlp_dc(x_k,Prob);
   end
   %[Prob.c_L Result.c_k Prob.c_U]
end
if nnLin > 0
   Result.Ax   = Prob.A*x_k;
   %[Prob.b_L Result.Ax Prob.b_U]
else
   Result.Ax = [];
end

% Lagrange multipliers
Result.v_k  = v_k;

Result.Iter     = iter_nlp;
%Result.FuncEv   = istat(4);
%Result.ConstrEv = istat(5);
%Result.GradEv   = istat(6);
%Result.HessEv   = istat(7);


Result.DUNDEE.ifail    = ifail;
Result.DUNDEE.lam      = v_k;
Result.DUNDEE.iter_nlp = iter_nlp;

optParam = Prob.optParam;

Result = StateDef(Result, x_k, Result.Ax, Result.c_k, optParam.xTol, ...
                  optParam.bTol, optParam.cTol, bl, bu, 0);

switch Result.Inform
   case 6
     ExitFlag=1;  % Too many iterations
   case 99999
     ExitFlag=2;  % Unbounded
   case {1,2}
     ExitFlag=4;  % Infeasible
   case {5,8}
     ExitFlag=3;  % Rank problem
   case {9,10,7,3,4}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.ExitFlag = ExitFlag;

switch(Result.Inform)
   
 case 0,
    if isempty(IntVars)
       ExitText = 'Optimal solution found';
    else
       ExitText = 'Optimal integer solution found';
    end
 case 1,
    ExitText = 'Root problem infeasible';
 case 2,
    ExitText = 'Integer infeasible';
 case 3,
    ExitText = 'Stack overflow - some i.f.s. obtained';
 case 4,
    ExitText = 'Stack overflow - no i.f.s. obtained';
 case 5,
    ExitText = 'SQP termination with rho < eps';
 case 6,
    ExitText = 'SQP termination with iter > max_iter';
 case 7,
    ExitText = 'Crash in user supplied routines';
 case 8,
    ExitText = 'Unexpected ifail from QP solvers';
 case 9,
    ExitText = 'Not enough REAL workspace or parameter error';
 case 10,
    ExitText = 'Not enough INTEGER workspace or parameter error';
 otherwise
    ExitText = 'Unknown ifail code from minlpBB';
end
Result.ExitText = ExitText;

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 021223 hkh  Use n = Prob.N. Default MaxIter 5000, stackmax 5000
% 021223 hkh  Use Prob.ConsPattern, no local ConsPattern
% 021223 hkh  Avoid calls to Prob.FUNCS.f etc. directly, could crash
% 021223 hkh  Compute pattern based on transpose of ConsPattern
% 021223 hkh  Errors in handling output c_k - corrected
% 021223 hkh  Change computation of Result.ExitFlag, Result.c_k, Result.Ax
% 030103 ango Change scaling parameters
% 030113 ango Add parameter DUNDEE.QPmin
% 030114 ango Add possibility to give multiplier vector DUNDEE.lam
% 030128 ango Revised comments
% 030129 ango Name changed MINLPbb->minlpBB, add safeguard: x_0
% 030129 ango Safer checking of bounds vectors
% 030206 hkh  Add optPar(20)=Nonlin. If 1, no special treatment of linear cons
% 030206 hkh  Add comments for optPar elements, set Result.x_0
% 030220 hkh  Must use mlp parameter in Mex, also maxf will later be needed
% 030220 hkh  Restrict mlp and kmax,  min(m, mlp), min(n,kmax)
% 030221 ango Added mlp&maxf parameters to MEX calls
% 030225 ango Change min value for mlp to m+1 or 3
% 031114 hkh  Change infty,BIG to 1E12 instead of 1E20, set default optPar(19)
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040209 ango Added parameter description: stackmax 
% 040330 hkh  optPar(9,10,11) was in wrong order compared to MEX code
% 040330 hkh  Set default stackmax to 10000
% 041202 hkh  Set default stackmax to 10000
% 041202 hkh  Revise calls to defblbu and StateDef, use different Order=0
% 041202 hkh  Use defblbu instead of messy code
% 041202 hkh  minlpbb fails if BIG < 1D20, correct this bug, safe guard BIG
% 041202 hkh  Code did not change inf to BIG, -inf to -BIG, done in defblbu
% 050914 med  iniSolve call changed
% 060714 med  Spelling updated
% 060814 med  FUNCS used for callbacks instead
% 060818 hkh  Use Prob.f_Low instead of -1E300 defining fLow, if QPmin not set 
% 060818 med  isnan checks removed for x_0
% 060829 med  LargeScale added in help
% 061210 med  MaxIter added to documentation
% 070222 hkh  Revised IntVars handling, use new format
% 070526 med  Output removed