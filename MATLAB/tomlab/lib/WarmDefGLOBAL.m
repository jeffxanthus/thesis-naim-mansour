% Initialization of structure Prob with the fields necessary for a warm
% start of the TOMLAB global optimization solvers.
%
% No checks are made, i.e. the routine will crash if the user is making
% the call with inappropriate input data i.e. fields that are undefined.
%
% function Prob = WarmDefGLOBAL(Solver,Prob,Result);
%
% INPUT:
%
% SolverName  Name of the TOMLAB GLOBAL optimization solver
%             Supported solvers:
%             FilMINT
%             multiMin
%             multiMINLP
%             glcDirect
%             glbDirect
%             glcCluster
%             ego
%             rbfSolve
%             arbfMIP
%             mipSolve
%             minlpSolve
%
% Prob        TOMLAB problem structure
% Result      The TOMLAB output structure from a previous run with
%             the same GLOBAL solver
%
% OUTPUT:
% Prob        Problem structure ready for warmstart

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 15, 2006.  Last modified July 28, 2011.

function Prob = WarmDefGLOBAL(SolverName,Prob,Result)

SolverName = upper(SolverName);

hasprob = isfield(Result,'Prob');

if hasprob
    if Result.Prob.N ~= Prob.N
        error('The number of decision variables are different');
    elseif Result.Prob.mLin ~= Prob.mLin
        error('The number of linear constraints are different');
    elseif Result.Prob.mNonLin ~= Prob.mNonLin
        error('The number of nonlinear constraints are different');
    end
end

switch lower(SolverName)
    case {'multimin'}
        Prob.WarmStart       = 1;
        Prob.multiMin        = Result.multiMin;
    case {'multiminlp'}
        Prob.WarmStart       = 1;
        Prob.multiMin        = Result.multiMin;
        Prob.multiMINLP      = Result.multiMINLP;
    case {'minlpsolve'}
        Prob.WarmStart       = 1;
        Prob.minlpSolve      = Result.minlpSolve;
    case {'filmint'}
        Prob.WarmStart       = 1;
        Prob.MINLP           = Result.MINLP;
    case {'mipsolve'}
        Prob.WarmStart       = 1;
        Prob.mipSolve        = Result.mipSolve;
    case {'glccluster'}
        Prob.WarmStart       = 1;
        Prob.MIP.fIP         = Result.f_k;
        Prob.MIP.xIP         = Result.x_k;
	Prob.Cluster         = Result.Cluster;
    case {'glb', 'glbdirect'}
        Prob.WarmStart = 1;
        Prob.glbDirect.WarmStartInfo = Result.glbDirect.WarmStartInfo;
    case {'glc', 'glcdirect'}
        Prob.WarmStart = 1;
        Prob.glcDirect.WarmStartInfo = Result.glcDirect.WarmStartInfo;
    case {'ego','arbfmip','rbfsolve','arbf'}
        Prob.WarmStart = 1;
        Prob.CGO.WarmStartInfo = Result.CGO.WarmStartInfo;
    otherwise
        error(sprintf('Unknown solver %s for warm start',SolverName));
end

% MODIFICATION LOG:
%
% 060715  hkh  Written
% 061004  ango Add ego
% 070514  med  Added size checks
% 080114  hkh  Added CGO solvers arbfMIP and rbfSolve
% 091002  hkh  Add CGO solver arbf (used in arbfMIP development) and minlpSolve
% 091002  hkh  Correcting check of glbDirect 
% 091004  hkh  Add solver mipSolve
% 091020  hkh  Added multiMINLP, revised multiMin
% 091118  hkh  Added FilMINT, use structure name MINLP
% 110615  hkh  Added glcCluster, set Prob.MIP.fIP = Result.f_k (fMin in solver)
% 110728  hkh  Added fields Result.Cluster to glcCluster warm start
