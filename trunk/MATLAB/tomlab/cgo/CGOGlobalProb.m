%function goProb = CGOGlobalProb(Prob, goName, glob_f, glob_g, glob_H,...
%   glob_c, glob_dc, glob_d2c, x_L, x_U, x_D, x_LL, x_UU, SCALE,...
%   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
%   backupSolver1, backupSolver2, PriSub)
%
% CGOGlobalProb generates a Prob structure for CGO subproblems
%
% Setting:
% goProb.optParam.MaxIter   = Prob.GO.MaxIter;
% goProb.optParam.MaxFunc   = Prob.GO.MaxFunc;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written April 10, 2008.    Last modified June 29, 2008.

function goProb = CGOGlobalProb(Prob, goName, glob_f, glob_g, glob_H,...
   glob_c, glob_dc, glob_d2c, x_L, x_U, x_D, x_LL, x_UU, SCALE,...
   dCon, dLin, IntVars, globalSolver, localSolver, bTol, cTol,...
   backupSolver1, backupSolver2, PriSub)

% HKH Note - x_U is never used

if Prob.simType == 1
   % Constraints are costly, avoid use in setting of global sub structure
   fprintf('Constraints are costly, avoid use in CGOGlobalProb\n');
   Prob.mNonLin   = 0;
   Prob.FUNCS.c   = [];
   Prob.FUNCS.dc  = [];
   Prob.FUNCS.d2c = [];
   Prob.c_L       = [];
   Prob.c_U       = [];
elseif Prob.simType == 2
   fprintf('Constraints both noncostly (used in CGOGlobalProb) and costly\n');
   % Constraints are both costly and noncostly. Noncostly used in expDesign
   % Noncostly constraints defined in Prob.CGO
   Prob.mNonLin     = Prob.CGO.mNonLin;
   Prob.ConsPattern = Prob.CGO.ConsPattern;
   Prob.FUNCS.c     = Prob.CGO.c;
   Prob.FUNCS.dc    = Prob.CGO.dc;
   Prob.FUNCS.d2c   = Prob.CGO.d2c;
   Prob.c_L         = Prob.CGO.c_L;
   Prob.c_U         = Prob.CGO.c_U;
end

dc = Prob.FUNCS.dc;

if SCALE > 0 & (dCon(1) > 0 | dLin > 0)
   if dLin > 0
      A   = Prob.A*diag(x_D);
      b_L = Prob.b_L - Prob.A*x_L;
      b_U = Prob.b_U - Prob.A*x_L;
   else
      A   = [];
      b_L = [];
      b_U = [];
   end
   if dCon(1) > 0
%function Prob = minlpAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
%                          IntVars, VarWeight, fIP, xIP, ...
%                          A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U,...
%                          x_min, x_max, f_opt, x_opt);
% function Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
%                            pSepFunc, fLowBnd, ...
%                            A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
%                            x_min, x_max, f_opt, x_opt);

      goProb             = minlpAssign(glob_f,glob_g,glob_H,[],x_LL,x_UU,goName,...
                           [],IntVars, [],[],[], ...
                           A, b_L, b_U, glob_c, glob_dc, glob_d2c,...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
      goProb.xD          = x_D;
      goProb.xL          = x_L;
      goProb.cNargin     = xnargin(Prob.FUNCS.c);
      if isempty(dc)
         goProb.dcNargin = 0;
      else
         goProb.dcNargin = xnargin(dc);
      end
      goProb.c           = Prob.FUNCS.c;
      goProb.dc          = Prob.FUNCS.dc;
      goProb.d2c         = Prob.FUNCS.d2c;
   else
      goProb             = minlpAssign(glob_f, glob_g,[],[],x_LL,x_UU,goName, ...
                           [], IntVars, [], [], [], A, b_L, b_U); 
   end

elseif dCon(1) == 0 & dLin == 0
   goProb                = minlpAssign(glob_f,glob_g,[],[],x_LL,x_UU,goName,[],IntVars); 
else
   goProb                = minlpAssign(glob_f, glob_g,[],[],x_LL,x_UU,goName, ...
                           [], IntVars, [], [], [], Prob.A, Prob.b_L, Prob.b_U, ...
                           Prob.FUNCS.c, Prob.FUNCS.dc, Prob.FUNCS.d2c, ...
                           Prob.ConsPattern, Prob.c_L, Prob.c_U);
end
goProb.SCALE              = SCALE;
goProb.SIGN               = 1;

d                         = length(x_D);
goProb.dDim               = d;
goProb.GO                 = Prob.GO;
goProb.GradTolg           = Prob.GradTolg;
goProb.GradTolH           = Prob.GradTolH;
goProb.GradTolJ           = Prob.GradTolJ;
goProb.MIP                = Prob.MIP;
goProb.uP                 = Prob.uP;
goProb.P                  = Prob.P;
goProb.Prilev             = PriSub;

% Solver special information
goProb.CONOPT             = Prob.CONOPT;
goProb.DUNDEE             = Prob.DUNDEE;
goProb.GENO               = Prob.GENO;
goProb.glbDirect          = Prob.glbDirect;
goProb.glcDirect          = Prob.glcDirect;
goProb.KNITRO             = Prob.KNITRO;
goProb.LGO                = Prob.LGO;
goProb.LSGRG2             = Prob.LSGRG2;
goProb.MISQP              = Prob.MISQP;
goProb.OQNLP              = Prob.OQNLP;
goProb.PENOPT             = Prob.PENOPT;
goProb.SOL                = Prob.SOL;
% Send user info
if isfield(Prob,'user')
   goProb.user            = Prob.user;
end
if isfield(Prob,'USER')
   goProb.USER            = Prob.USER;
end
% HKH Not needed:
%goProb.mLin               = dLin;
% HKH Needed:
goProb.mNonLin            = dCon(1);

% Define optParam parameters for globalSolver in goProb

solvType                  = checkType('minlp');
optParam                  = optParamDef(globalSolver,solvType,d,...
                                        dCon(1),dCon(1)+dLin);
goProb.optParam           = optParam;
% HKH good choice to use GOMaxFunc = Prob.GO.MaxFunc and
% and GOMaxIter = Prob.GO.MaxIter ???????????????????????
goProb.optParam.MaxIter   = Prob.GO.MaxIter;
goProb.optParam.MaxFunc   = Prob.GO.MaxFunc;
goProb.optParam.bTol      = bTol;
goProb.optParam.cTol      = cTol;
goProb.optParam.eps_Rank  = Prob.optParam.eps_Rank;
goProb.CGO.backupSolver1  = backupSolver1;
goProb.CGO.backupSolver2  = backupSolver2;

% Set LargeScale if user has set it in Prob
goProb.LargeScale         = Prob.LargeScale;

% Use goProb to define a local optimization structure ProbL
ProbL                     = goProb;

% Define optParam parameters for localSolver in ProbL
solvType                  = checkType('con');
optParam                  = optParamDef(localSolver,solvType,d,...
                                        dCon(1),dCon(1)+dLin);
ProbL.optParam            = optParam;
ProbL.optParam.MaxIter    = Prob.optParam.MaxIter;
ProbL.optParam.MaxFunc    = 2*Prob.optParam.MaxIter;
ProbL.optParam.bTol       = bTol;
ProbL.optParam.cTol       = cTol;
ProbL.optParam.eps_Rank   = Prob.optParam.eps_Rank;

% =========================================================================
% Take advantage of any derivatives in local search
% =========================================================================

if isempty(dc)
   cDiff                 = 6;
else
   cDiff                 = 0;
end
if isempty(glob_g)
   fDiff                 = 6;
else
   fDiff                 = 0;
end


switch lower(localSolver)
 case {'snopt','npsol','nlssol','minos'}
   if fDiff == 6
      if cDiff == 6
         DerLvl=0;
      else
         DerLvl=2;
      end
   elseif cDiff == 6
      DerLvl=1;
   else
      DerLvl=3;
   end
   ProbL.NumDiff        = fDiff;
   ProbL.ConsDiff       = cDiff;
   ProbL.SOL.optPar(39) = DerLvl;
   ProbL.SOL.optPar(12) = 1E-8; % Minor optimality tolerance
   %ProbL.SOL.optPar(10)= 1E-8;

 otherwise
   ProbL.NumDiff        = fDiff > 0;
   ProbL.ConsDiff       = cDiff > 0;
end

goProb.NumDiff        = fDiff > 0;
goProb.ConsDiff       = cDiff > 0;

% =========================================================================
% Save local structure ProbL in goProb.GO.ProbL for use for local search
% =========================================================================
goProb.GO.ProbL           = ProbL;


% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080410 hkh GOMaxFunc,GOMaxIter from Prob.GO instead of input variables
% 080416 hkh Set all solver fields into goProb from Prob
% 080416 hkh Also set goProb.MISQP from Prob.MISQP
% 080624 hkh Set LargeScale in sub structure using Prob.LargeScale
% 080629 hkh Handle SimType 1,2 with costly constraints Cc(x)
