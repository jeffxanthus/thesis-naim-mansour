% TOMLAB LSGRG2 NLP Solver
%
% function Result = tomRun('lsgrg2',Prob);
%
% LSGRG2 solves constrained nonlinear problems of the type:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, variable bounds
%    b_L <= A*x  <= b_U, linear constraints
%    c_L <= c(x) <= c_L, nonlinear constraints
%
%  where x_L, x_U, x are n*1-vectors, b_L,b_U are m1*1-vectors, A is a dense
%  or sparse m1*n matrix, and c_L, c_U, c(x) are m2*1-vectors.
%
% ------------------------------------------------------------------
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%    A          Linear constraints coefficient matrix.
%    x_L, x_U   Bounds on variables.
%    b_L, b_U   Bounds on linear constraints.
%    c_L, c_U   Bounds on nonlinear constraints.
%
%               For equality constraints (or fixed variables), set
%               e.g. b_L(k) == b_U(k).
%
%    PriLevOpt  Print level in optimizer and MEX interface. Set
%               Prob.LSGRG2.options.IPR to set optimizer print
%               level separately.
%
%    LargeScale Flag telling whether to treat the problem as sparse (1) or
%               dense. If set to 1, the user should also provide a sparse
%               0-1 matrix in Prob.ConsPattern giving the nonzero pattern.
%
%    MaxCPU     Maximum allowed time in seconds for the LSGRG2 run.
%               It is also possible to set this through the Prob.LSGRG2.options.MAXTIME
%               parameter, in which case Prob.MaxCPU is ignored. LSGRG2's
%               default value for MAXTIME is 1000 seconds.
%
% Prob.optParam Structure with optimization parameters. The following fields are used:
%
%    MaxIter    Maximum number of iterations. Default is 10000.
%
% Prob.LSGRG2   Structure with special fields for the LSGRG2 solver:
%
%    options    Structure array with options.
%               See the TOMLAB /OQNLP User's Guide for
%               instructions and examples.
%
%    PrintFile  Name of file to receive the LSGRG2 iteration and results
%               log. Independent of PriLevOpt.
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
%   bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%            Equality == 3;
%   cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%            Equality == 3;
%
%   v_k        Lagrangian multipliers (for bounds + dual solution vector).
%
%   ExitFlag   Exit status.
%
%   Inform     LSGRG2 information parameter.
%
%   rc               Reduced costs. If ninf=0, last m == -v_k.
%   Iter             Number of iterations.
%   FuncEv           Number of function evaluations.
%   GradEv           Number of gradient evaluations.
%   ConstrEv         Number of constraint evaluations.
%   QP.B             Basis vector in TOMLAB QP standard.
%   Solver           Name of the solver ('lsgrg2').
%   SolverAlgorithm  Description of the solver.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Dec 15, 2004.   Last modified Aug 14, 2006.

function Result = lsgrg2TL(Prob)

if nargin < 1
  error('lsgrg2TL needs the Prob structure as input');
end

Prob.solvType = 3; % NLP (CON) solver

Prob = iniSolve(Prob,3,1,1);

Result=ResultDef(Prob);
Result.Solver = 'LSGRG2';

LargeScale = DefPar(Prob,'LargeScale',0);

switch(LargeScale)
    case 0,
        Result.SolverAlgorithm = 'Dense LSGRG2';
    otherwise,
        Result.SolverAlgorithm = 'Sparse LSGRG2';
end

PriLev = DefPar(Prob,'PriLevOpt',0);

% LSGRG2 can probably have large big.
BIG=1E20;
[bl,bu,n,m1,m2] = defblbu(Prob,BIG,1);

% LSGRG2 might have a problem with "reversed" bounds 
idx = find(bl>bu);
if ~isempty(idx)
   % Any bounds where bl>bu are adjusted to equalities
   bl(idx)=bu(idx);
   if Prob.Warning > 0
      fprintf('\nWarning: Adjusting reversed bounds.\nTo disable this message, set Prob.Warning=0\n\n');
   end
end

m = m1+m2;

% Safe guarded starting point x_0:
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max(bl(1:n),min(bu(1:n),x_0)); 

% Function value at x_0
Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

% LSGRG2 substructure with solver-specific fields
LSGRG2 = DefPar(Prob,'LSGRG2',[]);
PrintFile = DefPar(LSGRG2, 'PrintFile', '');
options = DefPar(LSGRG2, 'options',[]);

MaxCPU = DefPar(Prob,'MaxCPU', 1000);
if ~isfield(options,'MAXTIME') & ~isempty(MaxCPU)
   options.MAXTIME = MaxCPU;
end

MaxIter = DefPar(Prob, 'MaxIter', 10000);
if ~isfield(options,'ITLIM') & ~isempty(MaxIter)
  options.ITLIM = MaxIter;
end

Prob.LSGRG2.options = options;

if LargeScale
   if issparse(Prob.A)
      A = Prob.A;
   else
      A = sparse(Prob.A);
   end

   if isempty(Prob.ConsPattern)
      ConsPattern = sparse( ones(m2,n) );
   else
      ConsPattern = sparse( Prob.ConsPattern );
   end
   
   nz = nnz(A)+nnz(ConsPattern)+n;
else
   A = full(Prob.A);
   ConsPattern = [];
   nz = (m+1)*n;
end

if m2 > 0
   % Determine the sparse problem structure
   if ~isempty(ConsPattern)
      [ix,iy]=find(ConsPattern);      
      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(ConsPattern),ix,iy);
   end
end

% The call
[Inform, x_k, f_k, c_k, v_k,inbind,redgr] = ...
   lsgrg2(n,m1,m2,...
   nz,ConsPattern,A,bl,bu,...
   x_0,PriLev,PrintFile,Prob);

Result.LSGRG2.ifail    = Inform;
Result.LSGRG2.rmults   = v_k;
Result.LSGRG2.inbind   = inbind;
Result.LSGRG2.redgr    = redgr;

% Recalculate final values - if not successful, some outputs from
% lsgrg2 may be wrong
Result.x_k = x_k;
Result.f_k = nlp_f(x_k,Prob);

if ~isempty(Prob.FUNCS.g), Result.g_k = nlp_g(x_k,Prob); else Result.g_k =[]; end
Result.H_k =[];

if m2 > 0
  Result.c_k = nlp_c(x_k,Prob);
  if ~isempty(Prob.FUNCS.dc)
    Result.cJac = nlp_dc(x_k,Prob);
  end
end

if m1 > 0
  Result.Ax = Prob.A*x_k;
else
  Result.Ax = [];
end

% Crude set of exit texts
ExitText = [];
ExitFlag = [];

switch(Inform)
 case {-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1}
  switch(Inform)
   case -17,   Result.LSGRG2.ExitText='PROBLEM_STRUCTURE';  
   case -16,   Result.LSGRG2.ExitText='ANAJAC_BAD_COL';
   case -15,   Result.LSGRG2.ExitText='BAD_USER_NNZ';
   case -14,   Result.LSGRG2.ExitText='ALLOCTBL_OVERFLOW';
   case -13,   Result.LSGRG2.ExitText='MISSING_GCOMPX';
   case -12,   Result.LSGRG2.ExitText='MISSING_PARSH';
   case -11,   Result.LSGRG2.ExitText='MISSING_GCOMP';
   case -10,   Result.LSGRG2.ExitText='BAD_OPTION_VALUE';
   case -9,    Result.LSGRG2.ExitText='BOUNDS_ERROR';
   case -8,    Result.LSGRG2.ExitText='LINEAR_VARS_ERROR';
   case -7,    Result.LSGRG2.ExitText='BAD_NOBJ';
   case -6,    Result.LSGRG2.ExitText='DIMENSION_ERROR';
   case -5,    Result.LSGRG2.ExitText='BAD_COMMAND';
   case -4,    Result.LSGRG2.ExitText='INTERNAL_ERROR';
   case -3,    Result.LSGRG2.ExitText='INVERT_FAILURE';
   case -2,    Result.LSGRG2.ExitText='INSFMEMORY';
   case -1,    Result.LSGRG2.ExitText='BADINPUT';
  end
  
  ExitText='Setup error';
  ExitFlag=10;
  
 case  0,
  Result.LSGRG2.ExitText='STATUS_NOT_SET';
  ExitText='Status not set';
  ExitFlag=11;
  
 case  1,
  Result.LSGRG2.ExitText='KTC';
  ExitText='Kuhn-Tucker conditions satisfied';
  ExitFlag = 0;
  
 case  2,    
  Result.LSGRG2.ExitText='FRACTCHG';
  ExitText='Fractional change in objective too small';
  ExitFlag = 0;
  
 case  3,    
  Result.LSGRG2.ExitText='ALLREMEDIES';
  ExitText='All remedies failed';
  ExitFlag=0; % ??
  
 case  4,
  Result.LSGRG2.ExitText='ITERATIONS';
  ExitText='Too many iterations';
  ExitFlag=1;
  
 case  5,
  Result.LSGRG2.ExitText='UNBOUNDED';
  ExitText = 'Problem is unbounded';
  ExitFlag = 2;
  
 case {6,7,8,9,10}
  switch(Inform)
   case  6,    Result.LSGRG2.ExitText='INFEASIBLE_KTC';
   case  7,    Result.LSGRG2.ExitText='INFEASIBLE_FRACTCHG';
   case  8,    Result.LSGRG2.ExitText='INFEASIBLE_ALLREMEDIES';
   case  9,    Result.LSGRG2.ExitText='INFEASIBLE_ITERATIONS';
   case 10,    Result.LSGRG2.ExitText='INFEASIBLE';
  end
  ExitText = 'Problem infeasible';
  ExitFlag = 4;
  
 case 11,    Result.LSGRG2.ExitText='SETUP_SUCCESS';
 case 12,    Result.LSGRG2.ExitText='SHUTDOWN_SUCCESS';
 case 13,    Result.LSGRG2.ExitText='REDOBJ_CONSTVIOL'; 
 case 14,    Result.LSGRG2.ExitText='REDGRA_NB_LE_0'; 
 case 15,    Result.LSGRG2.ExitText='XDOT_COLLEN'; 
 case 16,    Result.LSGRG2.ExitText='XSAXPY_COLLEN'; 
 case 17,    Result.LSGRG2.ExitText='GETBAS_INSFMEM'; 
 case 18,    Result.LSGRG2.ExitText='XPIVOT_COLLEN'; 
 case 19,    Result.LSGRG2.ExitText='CHUZQ_BADPIVOT'; 
 case 20,    Result.LSGRG2.ExitText='XPIVOT_BASIS_ILLCOND'; 
 case 21,    Result.LSGRG2.ExitText='XPIVOT_BASIS_SING'; 
 case 22,    Result.LSGRG2.ExitText='XPIVOT_INSFMEM'; 
 case 23,    Result.LSGRG2.ExitText='XPIVOT_OTHER_ERR'; 
 case 24,    Result.LSGRG2.ExitText='CONSBS_REINVERT_BC'; 
 case 25,    Result.LSGRG2.ExitText='CONSBS_BASIS_STRUCTURE'; 
 case 26,    Result.LSGRG2.ExitText='CONSBS_NOINVERT_SEARCH'; 
 case 27,    Result.LSGRG2.ExitText='CONSBS_BASIC_SLACK'; 
 case 28,    Result.LSGRG2.ExitText='CONSBS_JPIV_0'; 
 case 29,    Result.LSGRG2.ExitText='CONSBS_NB_TOOBIG'; 
 case 30,    Result.LSGRG2.ExitText='CONDNM_BAD_NBLOCKS'; 
 case 31,    Result.LSGRG2.ExitText='CONDNM_BAD_COLLEN'; 
 case 32,    Result.LSGRG2.ExitText='DIREC_UPDATE_ERR'; 
 case 33,    Result.LSGRG2.ExitText='PH0FAC_INVERT_FAILURE '; 
 case 34,    Result.LSGRG2.ExitText='PH0PIV_BAD_INDEX'; 
 case 35,    Result.LSGRG2.ExitText='PH0PIV_XPIVOT_FAILURE'; 
 case 36,    Result.LSGRG2.ExitText='PH0PIV_BAD_ICOLS1'; 
 case 37,    Result.LSGRG2.ExitText='PH0PIV_BAD_ICOLS2'; 
 case 38,    Result.LSGRG2.ExitText='OTHER_RUNTIME'; 
 
 case 39,
  Result.LSGRG2.ExitText='USER_TERMINATION'; 
  ExitText = 'Terminated by user function';
  ExitFlag = 12;
 
 case 40,
  Result.LSGRG2.ExitText='JAC_OVERFLOW'; 
 
 case 41,
  Result.LSGRG2.ExitText='OPTQUEST_ERROR'; 
  ExitText = 'OPTQUEST Error';
  ExitFlag = 13;
  
 case 42,
  Result.LSGRG2.ExitText='TIME_LIMIT_EXCEEDED';
  ExitText = 'Time limit exceeded';
  ExitFlag = 1;
  
otherwise,
  Result.LSGRG2.ExitText='';
  ExitText=['Unknown return code ' num2str(Inform) ];
  ExitFlag = -1;
end

if isempty(ExitText)
  Result.ExitText = deblank(Result.LSGRG2.ExitText);
else
  Result.ExitText = deblank(ExitText);
end

Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

% Multipliers
Result.v_k = v_k;

% Variable/constraint states
optParam = Prob.optParam;

Result = StateDef(Result, x_k, Result.Ax, Result.c_k, optParam.xTol, ...
                  optParam.bTol, optParam.cTol, bl, bu, 1);

Result = endSolve(Prob,Result);

% MODIFICATION LOG
% 
% 041220 frhe First version made
% 050110 med  Revised code and help
% 050112 med  OQNLP field changed to LSGRG2
% 050113 frhe Safed for empty MaxCPU and MaxIter
% 060814 med  FUNCS used for callbacks instead