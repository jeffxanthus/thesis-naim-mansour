% Initialization of structure Prob with the fields necessary for a warm
% start of a MINLP (BQPD, FILTERSQP) solver.
%
% No checks are made, i.e. the routine will crash if the user is making
% the call with inappropriate input data i.e. fields that are undefined.
%
% function Prob = WarmDefDUNDEE(Solver,Prob,Result);
%
% INPUT:
%
% SolverName  Solver name (FILTERSQP or BQPD)
% Prob        TOMLAB problem structure
% Result      The TOMLAB output structure from a previous run with
%             the same MINLP solver
%
% OUTPUT:
% Prob        Problem structure after change

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Jun 12, 2004.  Last modified May 14, 2007.

function Prob = WarmDefDUNDEE(SolverName,Prob,Result)

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

switch(upper(SolverName))
    case {'FILTERSQP','FILSQP','FILTER'}
        Prob.WarmStart = 1;
        Prob.DUNDEE.istat = Result.DUNDEE.istat;
        Prob.DUNDEE.lws   = Result.DUNDEE.lws;
        Prob.DUNDEE.lam   = Result.DUNDEE.lam;
        Prob.x_0          = Result.x_k;
        if hasprob
            if isfield(Result.Prob, 'ConsPattern')
                Prob.ConsPattern = Result.Prob.ConsPattern;
            end
            if isfield(Result.Prob, 'HessPattern')
                Prob.HessPattern = Result.Prob.HessPattern;
            end
        end
    case {'BQPD'}
        Prob.WarmStart = 1;
        Prob.DUNDEE.x     = Result.DUNDEE.x;
        Prob.DUNDEE.k     = Result.DUNDEE.k;
        Prob.DUNDEE.ls    = Result.DUNDEE.ls;
        Prob.DUNDEE.e     = Result.DUNDEE.e;
        Prob.DUNDEE.peq   = Result.DUNDEE.peq;
        Prob.DUNDEE.lp    = Result.DUNDEE.lp;
        Prob.x_0          = Result.x_k;
    otherwise
        error(sprintf('Unknown solver %s for warm start',SolverName));
end

% MODIFICATION LOG:
%
% 040607  med  Written
% 041013  med  x_0 added to Prob
% 041014  ango Some synonyms for filterSQP added
% 041202  med  Added ConsPattern and HessPattern
% 070514  med  Moved pattern and Prob check, added size check