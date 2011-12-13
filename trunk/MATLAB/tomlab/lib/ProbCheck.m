% Routine to check preset value in Prob and set all undefined values in Prob.
%
% function Prob = ProbCheck(Prob, Solver, solvType, probType);
%
% INPUT:
%  Prob     Problem structure
%  Solver   Solver name
%  solvType Solver type (or optType, if solver not known)
%  probType Problem type
%
% OUTPUT:
%  Prob     Problem structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Aug 5, 1999.    Last modified July 23, 2011.

function Prob = ProbCheck(Prob, Solver, solvType, probType)

if Prob.CHECK==1
   % ProbCheck has already checked this structure
   return
end

if nargin < 4
   probType=[];
   if nargin < 3
      solvType=[];
   end
end

if isempty(probType)
   probType=Prob.probType; 
end

if isempty(solvType) && isfield(Prob,'solvType')
   solvType=Prob.solvType;
end
% Guess that solvType the same as probType if not set
% Guess that solvType is constrained if not either probType or solvType set
if isempty(solvType)
    solvType=probType;
    if isempty(solvType)
        solvType=3;
    end
end
Prob.solvType=solvType;

% Guess that probType the same as solvType if not set
if isempty(probType)
    probType=solvType; 
end
Prob.probType=probType;
typeStr = checkType(probType);

% Check the struct and define any undefined items
if isempty(Prob.Solver.Name), Prob.Solver.Name=Solver; end

if isstruct(Prob.uP)
   fprintf('\nProb.uP must be [] or a numerical vector\n\n');
   error('ProbCheck: Prob.uP is a structure');
end

if strcmp('ode',typeStr)
   if ~isfield(Prob,'SolverODE'),     Prob.SolverODE=[]; end
   % Check ODE solution specific fields in Prob.ODE
   Prob = odeProbCheck(Prob);
   % Check ODE solver related fields
   if ~isfield(Prob.ODE,'LSODE'),     Prob.ODE.LSODE=[]; end
   if ~isfield(Prob.ODE,'RKSUITE'),   Prob.ODE.RKSUITE=[]; end
   if ~isfield(Prob.ODE,'ML'),        Prob.ODE.ML=[]; end
   % Check no fields related to parameter estimation, assume correct
end
if ~isfield(Prob,'Threads')
   Prob.Threads = 0;
elseif Prob.Threads > 1 
   if ~all(Prob.MATLAB >= [7,5])
      Prob.Threads = 0;
   end
end

if strcmp('exp',typeStr)
   if ~isfield(Prob.ExpFit,'p'),       Prob.ExpFit.p=[]; end
   if ~isfield(Prob.ExpFit,'wType'),   Prob.ExpFit.wType=[]; end
   if ~isfield(Prob.ExpFit,'eType'),   Prob.ExpFit.eType=double(1); end
   if ~isfield(Prob.ExpFit,'infCR'),   Prob.ExpFit.infCR=[]; end
   if ~isfield(Prob.ExpFit,'dType'),   Prob.ExpFit.dType=[]; end
   if ~isfield(Prob.ExpFit,'geoType'), Prob.ExpFit.geoType=[]; end
   if ~isfield(Prob.ExpFit,'qType'),   Prob.ExpFit.qType=[]; end
   if ~isfield(Prob.ExpFit,'sigType'), Prob.ExpFit.sigType=[]; end
   if ~isfield(Prob.ExpFit,'lambda'),  Prob.ExpFit.lambda=[]; end
   if ~isfield(Prob.ExpFit,'alpha'),   Prob.ExpFit.alpha=[]; end
   if ~isfield(Prob.ExpFit,'beta'),    Prob.ExpFit.beta=[]; end
   if ~isfield(Prob.ExpFit,'x0Type'),  Prob.ExpFit.xoType=[]; end
   if ~isfield(Prob.ExpFit,'sumType'), Prob.ExpFit.sumType=[]; end
end

mLin = Prob.mLin;
mNonLin = Prob.mNonLin;
N = Prob.N;
M = mLin + mNonLin;

% Cases LS, LLS, CLS, EXP, NTS
if any( strcmp('ls',typeStr)  | strcmp('lls',typeStr)  ...
      | strcmp('cls',typeStr) | strcmp('exp',typeStr)  ...
      | strcmp('nts',typeStr))
   if ~isfield(Prob.LS,'weightType'), Prob.LS.weightType=double(0); end
   if ~isfield(Prob.LS,'weightY'),    Prob.LS.weightY=[]; end
   if ~isfield(Prob.LS,'t'),          Prob.LS.t=[]; end
   if ~isfield(Prob.LS,'y'),          Prob.LS.y=[]; end
   if ~isfield(Prob.LS,'C'),          Prob.LS.C=[]; end
   if ~isfield(Prob.LS,'damp'),       Prob.LS.damp=[]; end
   if ~isfield(Prob.LS,'L'),          Prob.LS.L=[]; end
   if ~isfield(Prob.LS,'yUse'),       Prob.LS.yUse=1; end
   if ~isfield(Prob.LS,'SepAlg'),     Prob.LS.SepAlg=0; end

   % If empty FUNCS.f and FUNCS.fc, assume least squares problem
   if isempty(Prob.FUNCS.f) && isempty(Prob.FUNCS.fc)
      Prob.FUNCS.f = 'ls_f';
      Prob.FUNCS.g = 'ls_g';
      if strcmp('ls',typeStr)
         Prob.FUNCS.H = 'lls_H';
      else
         Prob.FUNCS.H = 'ls_H';
      end
   end
end

if ~isfield(Prob,'optParam')
   Prob.optParam = optParamDef(Solver,probType,N,N,M);
else
   % Check the optParam structure
   Prob.optParam = optParamSet(Prob.optParam,Solver,probType,N,N,M);
end

Prob.CHECK=1;

% MODIFICATION LOG:
%
% 990812 hkh  Expand QP fields.
% 990824 hkh  Expand QP fields with solver selections.
% 990909 hkh  Adding HessPattern, JacPattern and ConsPattern
% 000830 hkh  Fix for SOL subfields
% 000922 hkh  Correct the logic for global optimization, fields must be defined
% 001014 hkh  Added test on ConsDiff
% 010903 hkh  Now 63 parameters in optPar
% 011110 hkh  Add two new fields, GO and RBF, for global optimization
% 011112 hkh  Bug in definition of Prob.N, if a global opt problem defined as
%             being of another type
% 011213 hkh  Changed Prob.MIP.SC into SC, SI, semi-continuous and semi-integer
% 011226 hkh  Changed Prob.MIP.SOS1 and SOS2 to sos1 and sos2.
%             Initialize Prob.MIP.KNAPSACK = 0; not as empty
% 020103 hkh  Change field RBF to CGO
% 020409 hkh  Add field Mode and nState
% 020512 hkh  Check field DUNDEE for Fletcher and Leyffer Dundee solvers
% 020630 hkh  Check field PENSDP for Tomlab /PENSDP
% 021216 hkh  Add damp and L in field LS
% 030117 hkh  Change field PENSDP to PENOPT, for /PENSDP and /PENBMI
% 030117 hkh  Add check on field CheckNaN, default 0
% 030522 hkh  Add check on field FUNCS.fc and in check of least squares
% 030524 hkh  Add check on field FUNCS.gdc. Use double(0),double(1)
% 031129 hkh  Change fields AutoDiff to ADObj, ADCons
% 031201 hkh  Safeguard for ADObj <-1, ADCons <-1, ADMAT do not handle -2
% 040102 hkh  Compute fields mLin,mNonLin; Add checks on linear constraints
% 040115 hkh  Avoid working with mLin, mNonLin, N twice, set fields of missing
% 040115 hkh  Set fields mLin and mNonLin if missing
% 040303 hkh  Move definition of Prob.FUNCS.fc outside if-then-end-block
% 040425 hkh  Add field Prob.smallA + more, change order more like ProbDef
% 040506 hkh  conIx should be ConIx
% 040728 med  Keyboard removed
% 050302 hkh  Check if Prob.uP is a structure, if so, stop with error msg
% 050422 hkh  Check if Prob.fConstant is set, otherwise set 0
% 050503 hkh  Check if probType ODE, set missing fields
% 050605 hkh  Change to optParN = 65, for SNOPT 6
% 050606 hkh  Change to optParN = 71, for SNOPT 7
% 050616 hkh  Check on FUNCS.fg and FUNCS.cdc (and g_k)
% 050725 med  DUNDEE and PENOPT added optPar (options) vector default
% 050902 med  Added check on d2cPattern
% 051216 med  Check on FUNCS.rJ added
% 060814 med  FUNCS used for callbacks instead
% 080311 hkh  Set Prob.NumDiff, ConsDiff dependent on if user files are given
% 080416 hkh  Add field check for CONOPT, GENO, glbDirect, glcDirect, GP; 
% 080416 hkh  and KNITRO, LGO, LSGRG2, MILPSOLVE, OQNLP. Also MISQP.
% 080604 hkh  Not setting NumDiff, ConsDiff dep.on files, multilayer trouble
% 080606 med  Removed unnecessary field checks
% 080919 med  Replaced all checkType calls
% 090813 med  mlint check
% 110723 hkh  Check Prob.Threads defined. If incorrect Matlab version, set 0
