% Initialization of structure Prob with the fields necessary for a warm
% start of a SOL solver.
%
% No checks are made, i.e. the routine will crash if the user is making
% the call with inappropriate input data i.e. fields that are undefined.
%
% function Prob = WarmDefSOL(Solver,Prob,Result);
%
% INPUT:
%
% SolverName  Name of the SOL solver
% Prob        TOMLAB problem structure
% Result      The TOMLAB output structure from a previous run with
%             the same SOL solver
%
% OUTPUT:
% Prob        Problem structure after change

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written April 5, 2001.  Last modified May 14, 2007.

function Prob = WarmDefSOL(SolverName,Prob,Result)

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

switch SolverName
    case {'MINOS','SNOPT','SNOPT7','SQOPT','LP-MINOS','QP-MINOS'}
        Prob.WarmStart = 1;
        Prob.SOL.xs = Result.SOL.xs;
        Prob.SOL.hs = Result.SOL.hs;
        Prob.SOL.nS = Result.SOL.nS;
        if hasprob
            if isfield(Result.Prob, 'ConsPattern')
                Prob.ConsPattern = Result.Prob.ConsPattern;
            end
        end
    case {'NPSOL','NLSSOL'}
        Prob.WarmStart = 1;
        Prob.SOL.xs      = Result.SOL.xs;
        Prob.SOL.iState  = Result.SOL.iState;
        Prob.SOL.cLamda  = Result.SOL.cLamda;
        Prob.SOL.H       = Result.SOL.H;
        if hasprob
            if isfield(Result.Prob, 'ConsPattern')
                Prob.ConsPattern = Result.Prob.ConsPattern;
            end
        end
    case {'LPOPT','QPOPT','LSSOL'}
        Prob.WarmStart = 1;
        Prob.SOL.xs      = Result.SOL.xs;
        Prob.SOL.iState  = Result.SOL.iState;
end

% MODIFICATION LOG:
%
% 020208  hkh  Lambda changed to Lamda
% 041119  med  ConsPattern added
% 060228  ango Add SNOPT7
% 070514  med  Added size checks