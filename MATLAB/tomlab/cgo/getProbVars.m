% function [Prob ,MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, ... 
%           fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX, x_LL, x_UU, ...
%           x_DD, pDist, dLin, dCon, Percent, nSample, nTrial, AddMP,REPLACE,...
%           backupSolver1, backupSolver2, xOptS ] = ...
%           getProbVars(Prob, x_L,x_U,x_D, SCALE,globalSolver,IntVars, ...
%           Percent, nSample, nTrial, AddMP, REPLACE, GOMaxFunc, GOMaxIter, ...
%           maxFunc1, maxFunc2, maxFunc3, Solver)
%
% getProbVars picks up variables from the Prob structure, incl. Prob.optParam:
%   MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, 
%   fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX 
% 
% Default values (if []) are set for:
%   MaxIter, MaxFunc, IterPrint, epsX, PriSub
%
% dLin/dCon(1:2) are set using sizes of Prob.A and Prob.c_L/c_U
% dCon(1) Number of noncostly constraints, in Prob.CGO.c/c_L/c_U if
%         dCon(2) > 0, otherwise in Prob.FUNCS.c, Prob.c_L,c_U as usual
% dCon(2) Number of costly constraints, > 0 if Prob.simType > 0
%
% No adjustements of lengths in Prob.b_L/b_U/c_L/c_U (or Prob.CGO.c_L/c_U) 
% if lengths are less than dLin/dCon. cgoAssign/glcAssign takes care of this.
%
% getProbVars also defines values for variables from CGO structure:
%   Percent, nSample, AddMP, REPLACE
% Default values are dependent on dLin,dCon
%
% Default for Percent: (Really OK??????????? HKH)
%   Percent = -220 (ellipsoid 20%) if dLin+dCon(1)  > 0
%   Percent =  220 (ellipsoid 20%) if dLin+dCon(1) == 0
%
% Values are set for:
%   backupSolver1, backupSolver2, x_LL, x_UU, x_DD, pDist
%
% Scaled solution xOptS are defined if Prob.x_opt is given (optimal solution known)
%
% Prob.GO structure is set using the values given for:
% GOMaxFunc, GOMaxIter, maxFunc1, maxFunc2, maxFunc3
% The default values for these variables are dependent on the subsolvers used.
%


% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written April 10, 2008.    Last modified Nov 6, 2009.

function [Prob ,MaxCPU, PriLev, f_Low, MaxIter, MaxFunc, IterPrint, ... 
          fTol, fGoal, epsRank, bTol, cTol, PriSub, epsX, x_LL, x_UU, ...
          x_DD, pDist, dLin, dCon, Percent, nSample, nTrial, AddMP, REPLACE,...
          backupSolver1, backupSolver2, xOptS ] = ...
          getProbVars(Prob, x_L,x_U,x_D, SCALE,globalSolver,IntVars, ...
          Percent, nSample, nTrial, AddMP, REPLACE, GOMaxFunc, GOMaxIter, ...
          maxFunc1, maxFunc2, maxFunc3, Solver)

% Set possibly redefined SCALE into Prob
Prob.CGO.SCALE     = SCALE;

% Pick up input variables from Prob structure
MaxCPU    = Prob.MaxCPU;
PriLev    = Prob.PriLevOpt;          % Print level
f_Low     = Prob.f_Low;              % Lower bound on f
MaxIter   = Prob.optParam.MaxIter;   % Iterations used in global subopt
MaxFunc   = Prob.optParam.MaxFunc;   % Number of function evaluations
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration

fTol      = Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
fGoal     = Prob.optParam.fGoal;     % Goal f(x) for the optimization
epsRank   = Prob.optParam.eps_Rank;  % Rank tolerance in qr-decomposition
bTol      = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance
epsRank   = max(1E-14,epsRank);
if length(fGoal) > 1, fGoal = fGoal(1); end

% Tolerance used estimating components on bounds
epsX      = max(1E-3,fTol);
% epsX      = 1.7E-6; OLD value used in ARBFMIP

if isfield(Prob,'PriLevSub') 
   PriSub = Prob.PriLevSub;
else
   PriSub = [];
end
if isempty(PriSub), PriSub = 0; end

% Safeguard
if isempty(MaxFunc), MaxFunc = 300; end
if MaxFunc < 0
   MaxFunc = 300;
end
MaxFunc = min(MaxFunc,5000); % No safe guard needed, safety in mex

if isempty(MaxIter), MaxIter = 1000; end
if MaxIter < 0
   MaxIter = 1000;
end

if isempty(IterPrint), IterPrint = 1; end

% ----------------------------------------------------------------------------------------
d    = length(x_D);
dLin = size(Prob.A,1);
%if dLin > 0
%   Prob.b_L = [Prob.b_L;-inf*ones(dLin-length(Prob.b_L),1)];
%   Prob.b_U = [Prob.b_U; inf*ones(dLin-length(Prob.b_U),1)];
%end
% Nonlinear constraints may be both costly and noncostly
if Prob.simType == 0
   dCon(1) = Prob.mNonLin;
   % dCon(1) = max(length(Prob.c_L),length(Prob.c_U));
   dCon(2) = 0;
elseif Prob.simType == 1
   dCon(1) = 0;
   dCon(2) = Prob.mNonLin;
   % dCon(2) = max(length(Prob.c_L),length(Prob.c_U));
elseif Prob.simType == 2
   dCon(1) = Prob.CGO.mNonLin;
   dCon(2) = Prob.mNonLin;
else
   fprintf('Prob.simType %f\n',Prob.simType);
   error('getProbVars: Illegal value of Prob.simType');
end
% Default values for GO and CGO structure
if dLin + dCon(1) > 0
   if isempty(REPLACE),      REPLACE = 0; end
   % Default is constrained LH
   if isempty(Percent)      
      Percent = -6; 
      nSample = [];
      nTrial  = min(10000,d*500);
   end
   if isempty(IntVars)
      if isempty(maxFunc1),     maxFunc1  = 500*d; end
      if isempty(maxFunc2),     maxFunc2  = 1000*d; end
      if isempty(maxFunc3),     maxFunc3  = 1500*d; end
   else
      if isempty(maxFunc1),     maxFunc1  = 600*d; end
      if isempty(maxFunc2),     maxFunc2  = 1200*d; end
      if isempty(maxFunc3),     maxFunc3  = 2400*d; end
   end
   
   switch lower(globalSolver)
     case {'glccluster'}
        if isempty(IntVars)
           if isempty(GOMaxFunc),    GOMaxFunc = max(10000,2000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(5000,1000*d); end
        else
           if isempty(GOMaxFunc),    GOMaxFunc = max(20000,4000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(10000,2000*d); end
        end
     case {'glcdirect','glcfast','glcsolve'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(3000,300*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(3000,300*d); end
     case {'multimin'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
     case {'minlpbb'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
     otherwise
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
   end
else
   if isempty(REPLACE),      REPLACE = 5; end
   if isempty(Percent)
      if strcmpi(Solver,'ego')
         if d > 10
            Percent = 997;
         else
            Percent = 6;
            % if isempty(nSample), nSample = []; end
         end
      else
         % RBF solver should start with corner points, either 997 or 999
         %Percent = 999;
         Percent = 997;
      end
   end
   %if isempty(Percent),      Percent = 220; end     % Ellipsoid, 20%
   %if isempty(Percent),      Percent = -d; end
   %if isempty(Percent),      Percent = 997; end
   if isempty(IntVars)
      if isempty(maxFunc1),     maxFunc1  = 200+d*200; end
      if isempty(maxFunc2),     maxFunc2  = 0; end
      if isempty(maxFunc3),     maxFunc3  = maxFunc1; end
      %if isempty(maxFunc1),     maxFunc1  = 500*d; end
      %if isempty(maxFunc2),     maxFunc2  = 1000*d; end
      %if isempty(maxFunc3),     maxFunc3  = 1500*d; end
   else
      % RANDOM methods cannot be used, code does not handle constraints or Ints
      if Percent > 100 & Percent < 400, Percent = 997; end

      if isempty(maxFunc1),     maxFunc1  = 600*d; end
      if isempty(maxFunc2),     maxFunc2  = 1200*d; end
      if isempty(maxFunc3),     maxFunc3  = 2400*d; end
   end
   switch lower(globalSolver)
     case {'glccluster'}
        if isempty(IntVars)
           if isempty(GOMaxFunc),    GOMaxFunc = max(10000,2000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(5000,1000*d); end
        else
           if isempty(GOMaxFunc),    GOMaxFunc = max(20000,4000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(10000,2000*d); end
        end
     case {'glcdirect','glcfast','glcsolve'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(2000,200*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(2000,200*d); end
     case {'glbdirect','glbfast','glbsolve'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(2000,200*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(2000,200*d); end
     case {'multimin'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
     case {'minlpbb'}
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
     otherwise
        if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
        if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
   end
end

if isempty(AddMP)
   if any(Percent == [900,997,998,999])
      % Add the midpoint for the corner and adjacent corner strategies
      AddMP = 1;
   else
      AddMP = 0; 
   end
end

% Scaling parameters
if SCALE
   % General code for variable LOWER/UPPER
   % LOWER = 0;
   % UPPER = 1;
   % x_LL = LOWER*ones(d,1);
   % x_UU = UPPER*ones(d,1);
   % x_DD = x_UU-x_LL;
   % Efficient code 
   x_LL = zeros(d,1);
   x_UU = ones(d,1);
   x_DD = ones(d,1);
else
   x_LL = x_L;
   x_UU = x_U;
   x_DD = x_D;
end

pDist = norm(x_UU-x_LL);

%NHQ   This is taken care of in expDesign
% % % Initial Design Defaults
% % if Percent > 0 && Percent < 100
% %    if isempty(nSample)
% %       nSample = d+1;
% %    elseif nSample <= 0
% %       nSample = (d+1)*(d+2)/2;
% %    else
% %       nSample = max(d+1,nSample);
% %    end
% % elseif Percent == 0
% %    if isempty(nSample)
% %       %nSample = [];
% %    elseif nSample <= 0
% %       nSample = [];
% %    else
% %       nSample = max(d+1,nSample);
% %    end
% % elseif Percent < 0
% %    if isempty(nSample)
% %       %nSample = [];
% %    elseif nSample < 0
% %       nSample = (d+1)*(d+2)/2;
% %    elseif nSample == 0
% %       nSample = d+1;
% %    else
% %       nSample = max(d+1,nSample);
% %    end
% % end

if isfield(Prob.CGO,'backupSolver1')
   backupSolver1 = Prob.CGO.backupSolver1;
else
   backupSolver1 = [];
end
if isfield(Prob.CGO,'backupSolver2')
   backupSolver2 = Prob.CGO.backupSolver2;
else
   backupSolver2 = [];
end

if isempty(backupSolver1) 
% HKH USE backup solvers
if checkLicense('oqnlp')
   if isempty(IntVars)
      backupSolver1      = 'multiMin';
      backupSolver2      = 'oqnlp';
   elseif length(IntVars) == d
      backupSolver1      = 'oqnlp';
      backupSolver2      = 'glcDirect';
   else
      backupSolver1      = 'oqnlp';
      backupSolver2      = 'multiMin';
   end
else
   if length(IntVars) == d
      backupSolver1      = 'glcDirect';
      if checkLicense('minlpBB')
         backupSolver2      = 'minlpBB';
      else
         backupSolver2      = [];
      end
   else
      backupSolver1      = 'multiMin';
      backupSolver2      = 'glcDirect';
   end
end
end

% Use optimal value for printing, if given
if isempty(Prob.x_opt)
   x_opt = [];
else
   x_opt = Prob.x_opt;
   if any(size(x_opt)==1)
      x_opt = x_opt(:);
   else
      x_opt = x_opt(1,:)';
   end
   if length(x_opt)~=d, x_opt=[]; end
end
if SCALE
   if ~isempty(x_opt)
      xOptS   = (x_opt(:)-x_L)./x_D; 
   else
      xOptS   = [];
   end
else
   xOptS   = x_opt(:);
end

% Set values back in Prob.GO structure
Prob.GO.maxFunc1 = maxFunc1;
Prob.GO.maxFunc2 = maxFunc2;
Prob.GO.maxFunc3 = maxFunc3;
Prob.GO.MaxIter  = GOMaxIter;
Prob.GO.MaxFunc  = GOMaxFunc;

% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080414 hkh Revised comments. Is Ellipsoid strategy best as default???
% 080420 hkh Safe guard use of optimal solution in Prob.x_opt
% 080531 hkh Lower case solver names in Case statements, add MINLPBB
% 080626 hkh Use CGO.backupSolver1 and 2 if defined. Define default Percent
% 080629 hkh Define nCon(1:2) for noncostly and costly constraints
% 081105 hkh Add calling Solver last as input
% 081105 hkh Percent default set dependent on solver, corners algs for RBFs
% 081105 hkh Add input/output nTrial, set defaults for constrained problems
% 090426 hkh For IntVars, wrong default value of Percent
% 100222 hkh Use Percent=6 (LH) as default for box-bounded EGO
