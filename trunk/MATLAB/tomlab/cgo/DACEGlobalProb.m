%function DACEProb = DACEGlobalProb(Prob, daceName, dace_f, dace_g, dace_H,...
%   pEst, p0, pLow, pUpp, CGOLIB, globalSolver, localSolver, ...
%   backupSolver1, backupSolver2, PriSub)
%
% DACEGlobalProb generates a Prob structure for the DACE subproblem
%
% Setting:
% DACEProb.optParam.MaxIter   = Prob.GO.MaxIter;
% DACEProb.optParam.MaxFunc   = Prob.GO.MaxFunc;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 16, 2008.

% --------------------------------------------------------------
function DACEProb = DACEGlobalProb(Prob, daceName, dace_f, dace_g, dace_H,...
   pEst, p0, pLow, pUpp, CGOLIB, globalSolver, localSolver, ...
   backupSolver1, backupSolver2, PriSub)
% --------------------------------------------------------------

DebugPriLev = 0;  % PriLev in suboptimization, i.e. DACEProb.PriLevOpt

d = Prob.N;

if 1
   % Could use Fredriks assumption, correlation always between [0.01, 0.99]

   [xLL, xUU] = cgolib(208, CGOLIB.daceid);
   xLL        = log(xLL);
   xUU        = log(xUU);

else
   % or just set wide bounds for theta (log(alpha))
   % log(1E-10) =  -23  exp(10) = 22026
   % Dangerous to get stuck in local minima around -23
   xLL = -23*ones(d,1);
   xUU =  10*ones(d,1);
end

if pEst == 0
   x_L = xLL;
   x_U = xUU;
   x_0 = ones(d,1);
elseif pEst == 1
   x_L = [xLL;pLow];
   x_U = [xUU;pUpp];
   x_0 = [ones(d,1);p0(1)];
else
   x_L = [xLL;pLow];
   x_U = [xUU;pUpp];
   x_0 = [ones(d,1);p0];
end

DACEProb = conAssign(dace_f,dace_g,dace_H,[],x_L,x_U,daceName,x_0,[],[]);

% ---------------
% Used in dace_f:
% ---------------
DACEProb.d       = d;
DACEProb.pEst    = pEst;
if pEst == 0
   DACEProb.p    = p0*ones(d,1);
end
% ---------------------------------

DACEProb.GO                 = Prob.GO;
DACEProb.GradTolg           = Prob.GradTolg;
DACEProb.GradTolH           = Prob.GradTolH;
DACEProb.GradTolJ           = Prob.GradTolJ;
DACEProb.P                  = Prob.P;
DACEProb.Prilev             = PriSub;

% Solver special information
DACEProb.CONOPT             = Prob.CONOPT;
DACEProb.DUNDEE             = Prob.DUNDEE;
DACEProb.GENO               = Prob.GENO;
DACEProb.glbDirect          = Prob.glbDirect;
DACEProb.glcDirect          = Prob.glcDirect;
DACEProb.KNITRO             = Prob.KNITRO;
DACEProb.LGO                = Prob.LGO;
DACEProb.LSGRG2             = Prob.LSGRG2;
DACEProb.MISQP              = Prob.MISQP;
DACEProb.OQNLP              = Prob.OQNLP;
DACEProb.PENOPT             = Prob.PENOPT;
DACEProb.SOL                = Prob.SOL;
% MIP and user information not relevant for the DACE problem:
% Send user info
% if isfield(Prob,'user')
%    DACEProb.user            = Prob.user;
% end
% if isfield(Prob,'USER')
%    DACEProb.USER            = Prob.USER;
% end
%DACEProb.MIP                = Prob.MIP;
%DACEProb.uP                 = Prob.uP;

DACEProb.CGOLIB           = CGOLIB;
DACEProb.CGOLIB.TRANSFORM = 0;

% Use goProb to define a local optimization structure ProbL
ProbL                     = DACEProb;

DACEProb.optParam           = optParamDef(globalSolver,checkType('glb'),DACEProb.N,0,0);

% Define optParam parameters for globalSolver in DACEProb
DACEProb.PriLevOpt          = DebugPriLev;
DACEProb.optParam.IterPrint = DebugPriLev > 0;
DACEProb.optParam.MaxIter   = Prob.GO.MaxIter;
DACEProb.optParam.MaxFunc   = Prob.GO.MaxFunc;
DACEProb.optParam.eps_Rank  = Prob.optParam.eps_Rank;
DACEProb.CGO.backupSolver1  = backupSolver1;
DACEProb.CGO.backupSolver2  = backupSolver2;

% Define optParam parameters for localSolver in ProbL
solvType                  = checkType('con');
optParam                  = optParamDef(localSolver,solvType,d,0,0);
ProbL.optParam            = optParam;
ProbL.optParam.MaxIter    = Prob.optParam.MaxIter;
ProbL.optParam.MaxFunc    = 2*Prob.optParam.MaxIter;
ProbL.optParam.eps_Rank   = Prob.optParam.eps_Rank;

% =========================================================================
% Take advantage of any derivatives in local search
% =========================================================================

cDiff                    = 0;
if isempty(dace_g)
   fDiff                 = 6;
else
   fDiff                 = 0;
end

switch lower(localSolver)
 case {'snopt','npsol','nlssol','minos'}
   %if fDiff == 6
   %   DerLvl=2;
   %else
   %   DerLvl=3;
   %end
   ProbL.NumDiff        = fDiff;
   ProbL.ConsDiff       = cDiff;
   %ProbL.SOL.optPar(39) = DerLvl;
   ProbL.SOL.optPar(12) = 1E-8; % Minor optimality tolerance
   %ProbL.SOL.optPar(10)= 1E-8;

 otherwise
   ProbL.NumDiff        = fDiff > 0;
   ProbL.ConsDiff       = 0;
end

DACEProb.NumDiff        = fDiff > 0;
DACEProb.ConsDiff       = 0;

% =========================================================================
% Save local structure ProbL in DACEProb.GO.ProbL for use for local search
% =========================================================================
DACEProb.GO.ProbL           = ProbL;


% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080410 hkh GOMaxFunc,GOMaxIter from Prob.GO instead of input variables
% 080416 hkh Set all solver fields into DACEProb from Prob
% 080416 hkh Also set DACEProb.MISQP from Prob.MISQP
