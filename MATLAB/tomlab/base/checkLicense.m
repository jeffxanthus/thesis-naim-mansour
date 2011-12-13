% function [LicOK] = checkLicense(Solver)
%
% Solver    Name of TOMLAB solver (string)
% LicOK     Tells if 'Solver' is licensed for use or not (1/0)

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright 2007-2010 by Tomlab Optimization Inc., $Release: 7.5.0$
% Written Mar 2, 2007. Last modified Apr 12, 2010.

function [LicOK] = checkLicense(Solver)

if nargin < 1
    error('checkLicense requires one input argument')
end

[tomV,OS,TV] = tomlabVersion;

switch(lower(Solver))
    % 1. Base Module
    case {'glbdirect','glcdirect','lsei','milpsolve',...
            'qld','tfzero','tlsqr','tnnls','tomsol','base',...
            'multimin','multiminlp',...
            'consolve','nlpsolve','ucsolve','glccluster','slssolve',...
            'clssolve','glcfast','glbfast','glcsolve','glbscheckolve'}
        LicOK = TV(1);
        % 2. MINOS
    case {'minos','lpopt','qpopt'}
        LicOK = TV(2);
        % 3. NPSOL
    case {'npsol','lssol','nlssol'}
        LicOK = TV(3);
        % 4. SNOPT
    case {'snopt','sqopt','snopt7','sqopt7'}
        LicOK = TV(4);
        % 5. CGO
    case {'rbfsolve','ego','ego05','arbf','arbfmip','cgo'}
        LicOK = TV(5);
        % 6. PENSDP
    case {'pensdp'}
        LicOK = TV(6);
        % 7.   MINLP
    case {'minlpbb','bqpd','filsqp','filtersqp','miqpbb'}
        LicOK = TV(7);
        % 8.   Xpress
    case {'xpress-mp','xpress'}
        LicOK = TV(8);
        % 9.   CPLEX
    case {'cplex','cplexnet'}
        LicOK = TV(9);
        % 10.  PENBMI
    case {'penbmi'}
        LicOK = TV(10);
        % 11.  KNITRO
    case {'knitro'}
        LicOK = TV(11);
        % 12.  CONOPT
    case {'conopt'}
        LicOK = TV(12);
        % 13.  AMPL
    case {'ampl','amplfunc','spamfunc','amplqp'}
        LicOK = TV(13);
        % 14.  OQNLP
    case {'oqnlp'}
        LicOK = TV(14);
        % 15.  XA
    case {'xa'}
        LicOK = TV(15);
        % 16.  NLPQL
    case {'nlpql','nlpjob','dfnlp','nlpqlp'}
        LicOK = TV(16);
        % 17.  LGO
    case {'lgo'}
        LicOK = TV(17);
        % 18.  MAD
    case {'mad'}
        LicOK = TV(18);
        % 22.  MSNLP
    case {'msnlp','lsgrg2'}
        LicOK = TV(22);
        % 25.  GP
    case {'coplgp','gp'}
        LicOK = TV(25);
    case {'propt'}
        LicOK = TV(30);
    case {'gurobi'}
        LicOK = TV(32);
    otherwise
        LicOK = 0;
end

% MODIFICATION LOG:
%
% 070302 ango Wrote file
% 081006 med  Added multimin and cleaned up
% 090228 med  socs removed
% 090307 med  glcCluster added
% 090409 med  slsSolve and clsSolve added
% 100412 ango gurobi and propt added