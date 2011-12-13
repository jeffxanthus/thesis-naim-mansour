% tsprun:
%
% Run any of the predefined TSP problems from TSPLIB
%
% Presently 25 test problems are implemented.
%
% function tsprun(probNumber,SolverLP)
%
% INPUT
% probNumber:   Problem number
% SolverLP      The LP solver used by the solver salesman, default is lpSimplex
%               MINOS in TOMLAB /MINOS is recommended

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written June 15, 1999. Last modified Jul 31, 2006.

function tsprun(probNumber,SolverLP)

if nargin < 2
    SolverLP = [];
    if nargin < 1
        probNumber = 24;      % Problem number (=24 : ulysses16)
    end
end

Prob = probInit('tsp_prob', probNumber);
Name = Prob.Name;
m = size(Prob.TSP.C,1);
C = diag(NaN*ones(m,1)) + Prob.TSP.C;

fprintf([' Run problem : ' Name '\n'])

mu = zeros(m,1);         % Lagrangian multipliers

param.PriLev = 7;        % Print level
param.MaxIter = 100;     % Maximum # of iterations
param.DualGap = 0.001;   % Maximal duality gap in percent of dual objective
param.nPFS = 10;         % Compute PFS from iteration nPFS

param.ksi = 0.1;         % Step length parameter used by Polyak II
param.RedKsi = 3;        % Reduce ksi if dual fails to incr. in RedKsi iter.

[dual,pfs] = salesman(C,mu,param,SolverLP);

% MODIFICATION LOG:
%
% 990514 sto  First script version written by Susanne Timsjö
% 990704 hkh  Revised as a test function for TOMLAB v2.0
% 040816 hkh  Add extra input SolverLP, 4th parameter to salesman
% 060731 med  Modified to test prob format, new name tsprun