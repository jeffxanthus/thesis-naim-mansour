% cpxbiptest.m:
%
% Test of TOMLAB /CPLEX solving three binary integer linear 
% optimization problems calling CPLEX solver
%
% The test problem 1,2 have 1956 variables, 23 equalities and four inequalities
%
% Test problem 1, in bilp1.mat, is randomly generated, has several minima
% with optimal zero value
% Runs faster if avoiding the use of a cut strategy, and skipping presolve
%
% Test problem 2, in bilp2.mat, has unique minimum
% Runs faster if avoiding the use of presolve
%
% Test problem 3, in bilp1211.mat has 1656 variables, 23 equalities 
% and four inequalities. Runs very slow without the use of cuts.
% 
% function cpxbiptest(Cut, PreSolve, cpxControl)
%
% Cut       Value of the cut strategy control parameters, Default Cut = -1
%           -1 = Auto select of 1 or 2, 0 = No cuts, 1 = Conservative cuts, 
%           2  = Aggressive cut strategy
%
% PreSolve  Value of the presolve control parameter
%           0 = No presolve, 1 = presolve. Default = 1
%
%
% cpxControl The initial cpxControl structure. Here the user may set 
%            additional control parameters

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Sept 27, 2002. Last modified Feb 23, 2007.

function cpxbiptest(Cut, PreSolve, cpxControl)

if nargin < 3
   cpxControl = [];
   if nargin < 2
      PreSolve = [];
      if nargin < 1
         Cut = [];
end, end, end

if isempty(Cut),         Cut = -1; end
if isempty(PreSolve),    PreSolve = 1; end

% if Cut == 0, Cut = -1; end
% Cuts
cpxControl.CLIQUES    = Cut;
cpxControl.COVERS     = Cut;
cpxControl.DISJCUTS   = Cut;
cpxControl.FLOWCOVERS = Cut;
cpxControl.FLOWPATHS  = Cut;
cpxControl.FRACCUTS   = Cut;
cpxControl.GUBCOVERS  = Cut;
cpxControl.IMPLED     = Cut;
cpxControl.MIRCUTS    = Cut;

cpxControl.PREIND=PreSolve;

cpxControl

callback=[];

PriLev = 1;

Prob.P = 1;

% PROBLEM 1: COMBINATORIAL SEARCH
% -------------------------------

% Initial Problem, presumably many minima with zero value

load('bilp1.mat');

if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end
tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);

disp(' ')
disp('Press ENTER to continue')
pause
disp(' ')

% PROBLEM 2: Random Problem, presumably one minimum

load('bilp2.mat');

if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end
Prob.P = 2;

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);

% PROBLEM 3: 

load('bilp1211.mat');
if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end

disp(' ')
disp('Press ENTER to continue')
pause
disp(' ')
Prob.P = 3;

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);

% MODIFICATION LOG:
%
% 020927 hkh  Revision for CPLEX
% 021001 hkh  Add cuts and presolve
% 050113 med  Removed commented code
% 050117 med  mlint review
% 050209 med  Changed name from biptest to cpxbiptest
% 070223 hkh  Generate IntVars from noivars
