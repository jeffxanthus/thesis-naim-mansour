% Set values in structure optParam
%
% The fields in optParam store different optimization parameters
% that control the operation of the solvers, e.g. convergence tests
%
% function [optParam] = optParamSet(optParam,Solver,probType,nObj,nJac,m);
%
% INPUT:
%  optParam   Structure
%  Solver     The solver
%  probType   The type of problem/subproblem for the solver
%  nObj       Number of nonlinear objective variables
%  nJac       Number of nonlinear constraint variables
%  m          Number of rows in the constraint matrix (linear & nonlinear)
%
% OUTPUT:
%  optParam   Structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Sept 12, 1998. Last modified May 28, 2004.

function [optParam] = optParamSet(optParam,Solver,probType,nObj,nJac,m)

if isfield(optParam,'CHECK')
   if optParam.CHECK == 1, return; end
end

if nargin < 5
   m = [];
   if nargin < 4
      nJac = [];
      if nargin < 3
         nObj = [];
      end
   end
end

oP=optParamDef(Solver,probType, nObj, nJac, m);

if isempty(optParam)
   optParam = oP;
   optParam.CHECK = 1;
   return
end

if ~isfield(optParam,'PriLev')
   optParam.PriLev=oP.PriLev;
end
if ~isfield(optParam,'PriFreq')
   optParam.PriFreq=oP.PriFreq;
end
if ~isfield(optParam,'SummFreq')
   optParam.SummFreq=oP.SummFreq;
end
if ~isfield(optParam,'MinorPriLev')
   optParam.MinorPriLev=oP.MinorPriLev;
end
if ~isfield(optParam,'IterPrint')
   optParam.IterPrint=oP.IterPrint;
end
if ~isfield(optParam,'wait')
   optParam.wait=oP.wait;
end
if ~isfield(optParam,'MaxFunc')
   optParam.MaxFunc=oP.MaxFunc;
end
if ~isfield(optParam,'MaxIter')
   optParam.MaxIter=oP.MaxIter;
end
if ~isfield(optParam,'MajorIter')
   optParam.MajorIter=oP.MajorIter;
end
if ~isfield(optParam,'MinorIter')
   optParam.MinorIter=oP.MinorIter;
end
if ~isfield(optParam,'eps_f')
   optParam.eps_f=oP.eps_f;
end
if ~isfield(optParam,'eps_absf')
   optParam.eps_absf=oP.eps_absf;
end
if ~isfield(optParam,'eps_x')
   optParam.eps_x=oP.eps_x;
end
if ~isfield(optParam,'eps_dirg')
   optParam.eps_dirg=oP.eps_dirg;
end
if ~isfield(optParam,'eps_g')
   optParam.eps_g=oP.eps_g;
end
if ~isfield(optParam,'eps_Rank')
   optParam.eps_Rank=oP.eps_Rank;
end
if ~isfield(optParam,'EpsGlob')
   optParam.EpsGlob=oP.EpsGlob;
end
if ~isfield(optParam,'fGoal')
   optParam.fGoal=oP.fGoal;
end
if ~isfield(optParam,'fTol')
   optParam.fTol=oP.fTol;
end
if ~isfield(optParam,'xTol')
   optParam.xTol=oP.xTol;
end
if ~isfield(optParam,'bTol')
   optParam.bTol=oP.bTol;
end
if ~isfield(optParam,'cTol')
   optParam.cTol=oP.cTol;
end
if ~isfield(optParam,'MinorTolX')
   optParam.MinorTolX=oP.MinorTolX;
end
if ~isfield(optParam,'size_x')
   optParam.size_x=oP.size_x;
end
if ~isfield(optParam,'size_f')
   optParam.size_f=oP.size_f;
end
if ~isfield(optParam,'size_c')
   optParam.size_c=oP.size_c;
end
if ~isfield(optParam,'PreSolve')
   optParam.PreSolve=oP.PreSolve;
end
if ~isfield(optParam,'DerLevel')
   optParam.DerLevel=oP.DerLevel;
end
if ~isfield(optParam,'GradCheck')
   optParam.GradCheck=oP.GradCheck;
end
if ~isfield(optParam,'DiffInt')
   optParam.DiffInt=oP.DiffInt;
end
if ~isfield(optParam,'CentralDiff')
   optParam.CentralDiff=oP.CentralDiff;
end
if ~isfield(optParam,'QN_InitMatrix')
   optParam.QN_InitMatrix=oP.QN_InitMatrix;
end
if ~isfield(optParam,'splineSmooth')
   optParam.splineSmooth=oP.splineSmooth;
end
if ~isfield(optParam,'splineTol')
   optParam.splineTol=oP.splineTol;
end
if ~isfield(optParam,'BigStep')
   optParam.BigStep=oP.BigStep;
end
if ~isfield(optParam,'BigObj')
   optParam.BigObj=oP.BigObj;
end
optParam.CHECK = 1;

% MODIFICATION LOG
%
% 980918  hkh  Added line search parameters, incl. sigma and f_Low also.
% 980920  hkh  Added Penalty parameter for constrained problems
% 981005  hkh  Added new field items.
% 981009  hkh  Added new field LineSearch.MaxIter.
% 981026  hkh  Delete f_Low from optParam struct, put on Prob.f_Low
% 981027  hkh  Delete fields alg and subalg
% 981122  hkh  Add fields splineTol and splineSmooth
% 990910  hkh  Add parameter optParam.IterPrint
% 000709  hkh  Add parameter optParam.MinorIter
% 000911  hkh  Adding several parameters, ordering alphabetic
% 000925  hkh  Same order as structure, remove line search params
% 040102  hkh  First test if to skip further tests
% 040528  hkh  Set SolPAR(30) = 1000 for cplex, xpress-mp, make checks possible