% SOL solvers - Setting non-default parameters in optPar vector
%
% -----------------------------------------------------------------------
% function [optPar, SpecsFile, PrintFile, SummFile] = SOLSet(Solver, ...
%           probType, nObj, nJac, m, Prob)
% ------------------------------------------------------------------------
% INPUT:  
% Solver
% probType  String with type
% nObj      Number of nonlinear objective variables
% nJac      Number of nonlinear constraint variables
% m         Number of rows in the constraint matrix
%
% OUTPUT: 
% optPar

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Sep 3, 2000.     Last modified Oct 18, 2002.

function [optPar, SpecsFile, PrintFile, SummFile] = SOLSet(Solver, ...
          probType, nObj, nJac, m, Prob)

if nargin < 6
   error('Illegal call to SOLSet. Must have 6 input parameters')
end

if isempty(m),        m        = 70; end
if isempty(nJac),     nJac     = 70; end
if isempty(nObj),     nObj     = 70; end
if isempty(probType), probType = 3; end
if isempty(Solver),   Solver   = 'minos'; end

Solver = lower(Solver);

[optParDef] = SOLGet(Solver, probType, nObj, nJac, m);

optPar=Prob.SOL.optPar(:)';
optPar(length(optPar)+1:Prob.SOL.optParN)=-999;

SpecsFile=Prob.SOL.SpecsFile;
PrintFile=Prob.SOL.PrintFile;


if ~isempty(PrintFile)
   optPar(3)=9;
   %if optPar(1) <= 0, optPar(1)=11111; end
end
SummFile=Prob.SOL.SummFile;
if ~isempty(SummFile)
   optPar(4)=6;
   %if optPar(1) <= 0, optPar(1)=11111; end
end

% Set the correct parameters for the run

optParam  = Prob.optParam;
LineParam = Prob.LineParam;

if isempty(optParam),  return; end

if isempty(LineParam), LineParam.sigma=0.9; end
if isempty(LineParam.sigma) 
   if Solver(1)=='m' | Solver(1)=='M'
      LineParam.sigma=0.1;
   else
      LineParam.sigma=0.9;
   end
end

ixP = [9 10 11 12 22 27 30 35 36 41 5 6 42 43];

ix = find(optPar(ixP) == -999 & optParDef(ixP) > -100);

for i=1:length(ix)
    switch ix(i)
    case 1
      optPar(9)  = optParam.cTol;             % ROW Tolerance
    case 2
      optPar(10) = optParam.eps_x;            % Optimality Tol
    case 3
      optPar(11) = optParam.bTol;             % Feasibility Tol
    case 4
      optPar(12) = optParam.MinorTolX;        % Minor Optimality Tol
    case 5
      optPar(22) = LineParam.sigma; % Linesearch Tol
    case 6
      optPar(27) = optParam.eps_Rank;         % PIVOT TOL or Rank Tol
    case 7
      optPar(30) = optParam.MaxIter;          % Iterations
    case 8
      optPar(35) = optParam.MajorIter;        % Major Iterations
    case 9
      optPar(36) = optParam.MinorIter;        % Minor Iterations for QP/LP
                   % Feasibility Phase Iterations QP/LP QPOPT/LPOPT
    case 10
      optPar(41) = optParam.fTol;             % Function Precision
    case 11
      optPar(5)  = optParam.PriFreq;          % Print Frequency
    case 12
      optPar(6)  = optParam.SummFreq;         % Summary Frequency
      % optPar(39) = optParam.DerLevel;         % Derivative Level
    case 13
      optPar(42) = optParam.DiffInt;          % Difference Interval
    case 14
      optPar(43) = optParam.CentralDiff;      % Central Difference Interval
    end
end
iz = find(optPar(ixP) == optParDef(ixP));
if ~isempty(iz)
   optPar(ixP(iz))=-999;
end

% DEBUG: Print Level one
%optPar(1)=1;
%optPar(5)=1;
%optPar(8)=0;
% DEBUG: Use only function values
%optPar(39)=0;

if (strcmp(Solver(1:2),'sn') | Solver(1)=='m') & optPar(1)==-999
   % SNOPT or MINOS
   if optParam.IterPrint > 0
      optPar(1)=1;
   else
      PRI=optParam.PriLev;
      if PRI<=0
         %optPar(1)=0;
         % Make snopt and minos not creating any files
         optPar(1:3)=0;
      elseif PRI==1
         optPar(1)=1;
      elseif PRI==2
         optPar(1)=11111;
      else
         optPar(1)=111111;
      end
   end
elseif (Solver(1) == 'q' | Solver(1) == 'l') & optPar(1)==-999
   % QPOPT or LPOPT or LSSOL

   PRI=optParam.PriLev;
   if PRI<=0
      optPar(2)=0;
   elseif PRI==1
      optPar(2)=1;
   elseif PRI==2
      optPar(2)=5;
   elseif PRI==3
      optPar(2)=10;
   elseif PRI==4
      optPar(2)=20;
   else
      optPar(2)=30;
   end
end
if (Solver(2) == 'q' ) & optPar(2)==-999
   % SQOPT

   PRI=optParam.PriLev;
   if PRI<=0
      %optPar(2)=0;
      % Make sqopt not creating any files
      optPar(1:3)=0;
   elseif PRI==1
      optPar(2)=0;
   elseif PRI==2
      optPar(2)=1;
   else
      optPar(2)=10;
   end
end

%
% 000903 hkh  Written
% 010903 hkh  Avoid file creation for minos, snopt, sqopt if Print level 0
% 021214 hkh  Always make optPar row vector
% 041018 hkh  Allow LineParam.sigma [], set 0.1 if minos, else 0.1
