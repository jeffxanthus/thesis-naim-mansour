% TOMLAB OQNLP/MSNLP MINLP Solver
%
% Do not call this routine directly. Use oqnlpTL and msnlpTL instead.

% Fredrik Hellman, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Apr 20, 2005. Last modified Feb 23, 2007.

function Result = oqmsTL(Prob, Solver)

Prob.solvType = 12; % MINLP solver

Prob = iniSolve(Prob,12,1,1);

Result=ResultDef(Prob);

if Solver == 1
  Result.Solver = 'OQNLP';
else
  Result.Solver = 'MSNLP';
end

LargeScale = DefPar(Prob,'LargeScale',0);

switch(LargeScale)
 case 0,
  if Solver == 1
    Result.SolverAlgorithm = 'Dense Multistart GRG OQNLP 2.0';
  else
    Result.SolverAlgorithm = 'Dense Multistart GRG MSNLP 2.0';
  end
 otherwise,
  if Solver == 1
    Result.SolverAlgorithm = 'Sparse Multistart GRG OQNLP 2.0';
  else
    Result.SolverAlgorithm = 'Sparse Multistart GRG MSNLP 2.0';
  end
end

PriLev = DefPar(Prob,'PriLevOpt',0);

% OQNLP cannot have large BIG  
BIG=1E4;
[bl,bu,n,m1,m2] = defblbu(Prob,BIG,1);

% OQNLP might have a problem with "reversed" bounds 
idx = find(bl>bu);
if ~isempty(idx)
   % Any bounds where bl>bu are adjusted to equalities
   bl(idx)=bu(idx);
   if Prob.Warning > 0
      fprintf( ['\nWarning: Adjusting reversed bounds.\nTo disable this message, set Prob.Warning=0\n\n'] );
   end
end     
   
m = m1+m2;

% Safe guarded starting point x_0:
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max(bl(1:n),min(bu(1:n),x_0));

% Function value at x_0
Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

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
      error('oqmsTL: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

% OQNLP substructure with solver-specific fields
OQNLP = DefPar(Prob,'OQNLP',[]);

% iprint, not used anymore
%iprint    = DefPar(OQNLP,'iprint',[]);
%if isempty(iprint)
%    iprint = PriLev;
%end
iprint = [];

PrintFile = DefPar(OQNLP,'PrintFile','');

% StdFile, not used anymore. Should be removed from mex input argument list.
%StdFile   = DefPar(OQNLP,'StdFile','oqnlpstd.txt');
StdFile = '';

options = DefPar(OQNLP,'options',[]);
if ~isfield(options,'ITERATION_LIMIT')
   options.ITERATION_LIMIT = DefPar(Prob.optParam,'MaxIter',[]);
end

MaxCPU = DefPar(Prob,'MaxCPU',inf);
if ~isfield(options,'MAXTIME') & isfinite(MaxCPU)
   options.MAXTIME = MaxCPU;
end

Prob.OQNLP.options = options;

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

if ~isempty(IntVars)
    if Prob.simType > 0
        if isempty(Prob.FUNCS.gdc) | (Prob.ConsDiff > 0 ) | (Prob.NumDiff > 0)
            Prob.FUNCS.gdc = 'oqnlp_gdc';  % Special interface routine
            Prob.MIP.IntVars = IntVars;   % To avoid doing tests on IntVars
            Prob.CheckNaN = 1;
            Prob.ConsDiff = 0;
            Prob.NumDiff = 0;
            global NARG
            NARG(11) = 2;
            if isempty(ConsPattern)
                Prob.ConsPattern = ones(m2,n);
                Prob.ConsPattern(:,IntVars) = 0;
            else
                Prob.ConsPattern(:,IntVars) = 0;
            end
        end
    else
        if isempty(Prob.FUNCS.g) | (Prob.NumDiff > 0)
            Prob.FUNCS.g = 'oqnlp_g';      % Interface routine for gradient
            Prob.MIP.IntVars = IntVars;   % To avoid doing tests on IntVars
            Prob.CheckNaN = 1;
            global NARG
            NARG(2) = 2;
        end
        if isempty(Prob.FUNCS.dc) | (Prob.ConsDiff > 0 )
            Prob.FUNCS.gdc = 'oqnlp_dc';   % Interface routine for Jacobian
            Prob.MIP.IntVars = IntVars;   % To avoid doing tests on IntVars
            Prob.CheckNaN = 1;
            global NARG
            NARG(5) = 2;
            if isempty(ConsPattern)
                Prob.ConsPattern = ones(m2,n);
                Prob.ConsPattern(:,IntVars) = 0;
            else
                Prob.ConsPattern(:,IntVars) = 0;
            end
        end
    end
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
if Solver == 1
  [rc, x_k, f_k, c_k, v_k,inbind,redgr] = ...
      oqnlp(n,m1,m2,...
            nz,ConsPattern,A,bl,bu,...
            x_0,double(IV), PriLev,PrintFile,StdFile,iprint,Prob);
else
  [rc, x_k, f_k, c_k, v_k,inbind,redgr] = ...
      msnlp(n,m1,m2,...
            nz,ConsPattern,A,bl,bu,...
            x_0,double(IV), PriLev,PrintFile,StdFile,iprint,Prob);
end

% New type of return code, needs some bitwise arithmetics
if rc >= 0
   rc_loc = bitand(rc,2^16-1); % Filter out the low bits, which is the local solution status
   rc_alg = bitshift(rc,-16);
else
   rc_loc = -bitand(-rc,2^16-1);
   rc_alg = bitshift(-rc,-16);
end

%       i = oqnlpLocalSolStatus;
%       if(i < 0) {
%          i = -i;
%          info = (oqnlpAlgTermination << 16) + i;
%          info = -info;
%       }
%       else info = (oqnlpAlgTermination << 16) + i;

Inform   = rc_loc;
% "Raw" outputs from oqnlp:
Result.OQNLP.ifail    = rc;
Result.OQNLP.rc_loc   = rc_loc;
Result.OQNLP.rc_alg   = rc_alg;
Result.OQNLP.rmults   = v_k;
Result.OQNLP.inbind   = inbind;
Result.OQNLP.redgr    = redgr;

% Recalculate final values - if not successful, some outputs from
% oqnlp may be wrong
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

%
% Two status codes: Algorithm status and local solution status. The
% algorithm status tells us why the multi start algorithm stopped
% searching. The local solution status gives us information on the
% status of the best local solution. These two are combined into
% one big exit text.
%

ExitFlag = [];
AlgExitText = [];
AlgShortexitText = [];
LocExitText = [];

switch(rc_alg)
 case 1,
  AlgExitText = 'Iteration limit exceeded.';
 case 2,
  AlgExitText = 'Time limit exceeded.';
 case 3,
  AlgExitText = 'Locals limit exceeded.';
 case 4,
  AlgExitText = 'Solver call exceeded.';
 case 5,
  AlgExitText = 'Local solution improvement criterion met.';
 case 6,
  AlgExitText = 'User termination.';
 case 7,
  AlgExitText = 'Discrete optimum found.';
end

switch(rc_loc)
   case {-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1}
      switch(Inform)
         case -17,   LocShortExitText='PROBLEM_STRUCTURE';  
         case -16,   LocShortExitText='ANAJAC_BAD_COL';
         case -15,   LocShortExitText='BAD_USER_NNZ';
         case -14,   LocShortExitText='ALLOCTBL_OVERFLOW';
         case -13,   LocShortExitText='MISSING_GCOMPX';
         case -12,   LocShortExitText='MISSING_PARSH';
         case -11,   LocShortExitText='MISSING_GCOMP';
         case -10,   LocShortExitText='BAD_OPTION_VALUE';
         case -9,    LocShortExitText='BOUNDS_ERROR';
         case -8,    LocShortExitText='LINEAR_VARS_ERROR';
         case -7,    LocShortExitText='BAD_NOBJ';
         case -6,    LocShortExitText='DIMENSION_ERROR';
         case -5,    LocShortExitText='BAD_COMMAND';
         case -4,    LocShortExitText='INTERNAL_ERROR';
         case -3,    LocShortExitText='INVERT_FAILURE';
         case -2,    LocShortExitText='INSFMEMORY';
         case -1,    LocShortExitText='BADINPUT';
      end
      
    LocExitText='Setup error.';
    ExitFlag=10;
    
 case  0,    LocShortExitText='STATUS_NOT_SET';
    LocExitText='Status not set.';
    ExitFlag=11;
    
 case  1,
    LocShortExitText='KTC';
    LocExitText='Optimal solution found.';
    ExitFlag = 0;
    
 case  2,    LocShortExitText='FRACTCHG';
  LocExitText='Fractional change in objective too small.';
  ExitFlag = 0;
  
 case  3,    LocShortExitText='ALLREMEDIES';
  LocExitText='All remedies failed.';
  ExitFlag=0; % ??
  
 case  4,
  LocShortExitText='ITERATIONS';
  LocExitText='Too many iterations.';
  ExitFlag=1;
  
 case  5,
  LocShortExitText='UNBOUNDED';
  LocExitText = 'Problem is unbounded.';
  ExitFlag = 2;
  
    
 case {6,7,8,9,10}
    switch(Inform)
       case  6,    LocShortExitText='INFEASIBLE_KTC';
       case  7,    LocShortExitText='INFEASIBLE_FRACTCHG';
       case  8,    LocShortExitText='INFEASIBLE_ALLREMEDIES';
       case  9,    LocShortExitText='INFEASIBLE_ITERATIONS';
       case 10,    LocShortExitText='INFEASIBLE';
    end
    LocExitText = 'Problem infeasible.';
    ExitFlag = 4;
    
 case 11,    LocShortExitText='SETUP_SUCCESS';
 case 12,    LocShortExitText='SHUTDOWN_SUCCESS';
 case 13,    LocShortExitText='REDOBJ_CONSTVIOL'; 
 case 14,    LocShortExitText='REDGRA_NB_LE_0'; 
 case 15,    LocShortExitText='XDOT_COLLEN'; 
 case 16,    LocShortExitText='XSAXPY_COLLEN'; 
 case 17,    LocShortExitText='GETBAS_INSFMEM'; 
 case 18,    LocShortExitText='XPIVOT_COLLEN'; 
 case 19,    LocShortExitText='CHUZQ_BADPIVOT'; 
 case 20,    LocShortExitText='XPIVOT_BASIS_ILLCOND'; 
 case 21,    LocShortExitText='XPIVOT_BASIS_SING'; 
 case 22,    LocShortExitText='XPIVOT_INSFMEM'; 
 case 23,    LocShortExitText='XPIVOT_OTHER_ERR'; 
 case 24,    LocShortExitText='CONSBS_REINVERT_BC'; 
 case 25,    LocShortExitText='CONSBS_BASIS_STRUCTURE'; 
 case 26,    LocShortExitText='CONSBS_NOINVERT_SEARCH'; 
 case 27,    LocShortExitText='CONSBS_BASIC_SLACK'; 
 case 28,    LocShortExitText='CONSBS_JPIV_0'; 
 case 29,    LocShortExitText='CONSBS_NB_TOOBIG'; 
 case 30,    LocShortExitText='CONDNM_BAD_NBLOCKS'; 
 case 31,    LocShortExitText='CONDNM_BAD_COLLEN'; 
 case 32,    LocShortExitText='DIREC_UPDATE_ERR'; 
 case 33,    LocShortExitText='PH0FAC_INVERT_FAILURE '; 
 case 34,    LocShortExitText='PH0PIV_BAD_INDEX'; 
 case 35,    LocShortExitText='PH0PIV_XPIVOT_FAILURE'; 
 case 36,    LocShortExitText='PH0PIV_BAD_ICOLS1'; 
 case 37,    LocShortExitText='PH0PIV_BAD_ICOLS2'; 
 case 38,    LocShortExitText='OTHER_RUNTIME'; 
 
 case 39,
  LocShortExitText='USER_TERMINATION'; 
  LocExitText = 'Terminated by user function.';
  ExitFlag = 12;
 
 case 40,    LocShortExitText='JAC_OVERFLOW'; 
 
 case 41,
  LocShortExitText='OPTQUEST_ERROR'; 
  LocExitText = 'OPTQUEST Error.';
  ExitFlag = 13;
  
 case 42,
  LocShortExitText='TIME_LIMIT_EXCEEDED';
  LocExitText = 'Time limit exceeded.';
  ExitFlag = 1;

 case 43,
  LocShortExitText = 'FOUND_FEASIBLE';
  LocExitText = 'Feasible solution found.';
  ExitFlag = 0;
  
 otherwise,
  LocShortExitText='';
  LocExitText=['Unknown return code ' num2str(Inform) ' (' num2str(rc) ').' ];
  ExitFlag = -1;
end


if isempty(LocExitText)
  LocExitText = deblank(LocShortExitText);
else
  LocExitText = deblank(LocExitText);
end

if isempty(AlgExitText)
  AlgExitText = 'Unknown.';
end

Result.OQNLP.AlgExitText = AlgExitText;
Result.OQNLP.LocExitText = LocExitText;

Result.ExitText = sprintf('%s %s\n%s %s', ...
                          'Sol. status:', ...
                          LocExitText, ...
                          'Alg. status:', ...
                          AlgExitText);

Result.ExitFlag = ExitFlag;
Result.Inform   = rc_loc;

% Constraint values
% if m1>0, Result.Ax  = c_k(1:m1); else Result.Ax=[]; end
% Result.c_k = c_k(m1+1:m);
% 
% Result.cJac = nlp_dc(x_k,Prob);

% Multipliers
Result.v_k = v_k;

% Variable/constraint states
optParam = Prob.optParam;

Result = StateDef(Result, x_k, Result.Ax, Result.c_k, optParam.xTol, ...
                  optParam.bTol, optParam.cTol, bl, bu, 1);

Result=endSolve(Prob,Result);

% MODIFICATION LOG
% 
% 030414 ango Wrote file
% 030904 ango x/b/cStates now correctly returned
% 031118 ango LargeScale handling slightly changed
% 031210 ango Change comment on calling syntax
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040128 hkh  Safeguard for patological equalities in A, rhs > abs(BIG)
% 040507 ango Fixed for bit-combined return codes from OQNLP solver
% 040526 hkh  Wrong call to iniSolve, only 1st order information needed
% 040602 ango Added Prob.MaxCPU handling. 
% 040928 ango Safeguard for reversed bounds.
% 041112 hkh  Use oqnlp_gdc for automatic avoidance of IntVars for num.diff
% 041112 frhe Added help text about PrintFile.
% 041202 hkh  Revise calls to defblbu and StateDef, avoid vector reshuffling
% 041202 hkh  Unnecessary tests removed
% 050112 frhe Status codes handled totally different now.
% 050122 hkh  Inform used, but not set as rc_loc
% 050420 frhe New file: oqmsTL.m. Based on oqnlpTL.m. To be used by
%             both msnlp and oqnlp.
% 060814 med  FUNCS used for callbacks instead
% 070222 hkh  Revise IntVars handling, use new format
% 070223 hkh  Boolean input to mex must be given as double vector
