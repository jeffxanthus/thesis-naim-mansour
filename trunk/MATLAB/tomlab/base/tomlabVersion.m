% function [tomV,OS] = tomlabVersion
%
% tomV Tomlab Version
% OS   Operating system
%
% Returns the current Tomlab version

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 1, 2000.     Last modified Nov 30, 2011.

function [tomV,OS,TV] = tomlabVersion

% Tomlab v3.1: Version any of MINI, v3.0, SOL, NPSOL, SNOPT, CGO
%
% Tomlab v3.2 and v4.0:
% Version any of Base Module, MINOS, NPSOL, SNOPT, SOL, CGO, PENSDP, MINLP,
%                Xpress, CPLEX, PENBMI
%
% Values =1 in TV vector if license for module is OK. Indices described below:
%
% TV   Tomlab module
% 1.   Base Module
% 2.   MINOS
% 3.   NPSOL
% 4.   SNOPT
% 3+4. SOL
% 5.   CGO
% 6.   PENSDP
% 7.   MINLP
% 8.   Xpress
% 9.   CPLEX
% 10.  PENBMI
% 11.  KNITRO
% 12.  CONOPT
% 13.  AMPL
% 14.  OQNLP
% 15.  XA
% 16.  NLPQL
% 17.  LGO
% 18.  MAD
% 22.  MSNLP
% 25.  GP

[x1,x2,x3,x4,x5,x6,tomV,OS]=tomlablic(1);

NoOfProducts = 40;

if ischar(tomV)
    TV    = [1; zeros(NoOfProducts,1)];
    TV(1) = 1;
    if findstr('v3',tomV)
        TV(2) = 1;
    elseif findstr('SO',tomV)
        TV(2:4) = 1;
    elseif findstr('SN',tomV)
        TV(2) = 1;
        TV(4) = 1;
    elseif findstr('NP',tomV)
        TV(2) = 1;
        TV(3) = 1;
    end
    if findstr('CGO',tomV)
        TV(5) = 1;
    end
else
    TV   = [tomV(2:end);zeros(NoOfProducts-length(tomV),1)];
    tomV = tomlablic(2);
end

% Filter out components unavailable on current platform
switch(computer)
    %case 'PCWIN'  % WIN32, all available
    %   return;
    case 'PCWIN64' % WIN64
        TV( [8,15,20,21,23,24] ) = 0;
    case 'GLNX86'  % LNX32
        TV( [15,21,23,24] ) = 0;
    case 'GLNXA64' % LNX64
        TV( [8,14,15,20,21,23,24] ) = 0;
    case 'MAC' % MAC OSX PPC
        TV( [6,8,9,10,11,12,14,15,20,21,23,24] ) = 0;
    case 'MACI' % MAC OSX Intel
        TV( [6,8,10,12,14,15,20,21,23,24] ) = 0;
    case 'MACI64' % MAC OS X Intel x86_64
        TV( [6,8,10,14,15,20,21,23,24] ) = 0;
end

% MODIFICATION LOG:
%
% 070306 ango Add KNITRO for Win64
% 070307 ango Fix Win64 Boeing wrong index, split MAC, MACI
% 070410 ango Add CONOPT, LNX64
% 070911 ango Add CONOPT, WIN64
% 080108 ango Add CPLEX for MACI
% 090228 med  socs removed
% 110504 ango Win64 PENOPT added
% 111130 ango Win64 OQNLP added 
