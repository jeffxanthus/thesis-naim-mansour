% Test of TOMLAB /CPLEX.
%
% Running a generalized assignment problem from Wolsey 1998, 
% Section 9.6, pp159.
%
% Define the linear sos1 constraints explicitly
%
% function x = cpxTest3(DEFPARAM)
%
% DEFPARAM  If =1 use default parameters, presolve, cuts, Dual simplex, 
%           If =0 no presolve or cuts, Primal simplex (default)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Oct 28, 1999.   Last modified Feb 23, 2007.

function x = cpxTest3(DEFPARAM)

if nargin < 1
   DEFPARAM = [];
end

if isempty(DEFPARAM), DEFPARAM = 0; end
  
% Global Tomlab variables that determine number of elements printed
global MAX_x MAX_c


A=[95  1 21 66 59; ...
   54 53 44 26 60; ...
    3 91 43 42  5; ...
   72 30 56 72  9; ...
   44  1 71 13 27; ...
   20 99 87 52 85; ...
   72 96 97 73 49; ...
   75 82 83 44 59; ...
   68  8 87 74  4; ...
   69 83 98 88 45];

C=-[110  16 25 78 59; ...
     65  69 54 28 71; ...
     19  93 45 45  9; ...
     89  31 72 83 20; ...
     62  17 77 18 39; ...
     37 115 87 59 97; ...
     89 102 98 74 61; ...
     78  96 87 55 77; ...
     74  27 99 91  5; ...
     88  97 99 99 51];

b = [91 87 109 88 64]';

MIP=1;

[c,x_L,x_U,b_L,b_U,A]=abc2gap(A,b,C,0);

IntVars=[1:length(c)];
MAX_x = length(c);

% PriLev > 0 gives printing in this file
% Increasing PriLev gives more output

PriLev = 1;

PI      = [];
SC      = [];
SI      = [];
sos1    = [];
sos2    = [];

Prob.P = 3;  % Initialize Prob with problem number 3

callback=[]; % Use default callbacks

if PriLev > 0
   PriLev
end

if DEFPARAM
   % Empty cpxControl (default values) will give fastest execution
   cpxControl = [];
else
   % Setting Cplex Parameters (will slow down the computations)
   % Change to Primal simplex for root and subsolutions
   % 0 default 1 = Primal 2 = Dual (default) 3 = Network 4 = Barrier
   cpxControl.NODEALG    = 1;
   cpxControl.SUBALG     = 1;
   % No PreSolve
   cpxControl.PREIND     = 0;
   cpxControl.AGGIND     = 0;
   % Disable the cuts
   cpxControl.CLIQUES    = -1;
   cpxControl.COVERS     = -1;
   cpxControl.DISJCUTS   = -1;
   cpxControl.FLOWCOVERS = -1;
   cpxControl.FLOWPATHS  = -1;
   cpxControl.FRACCUTS   = -1;
   cpxControl.GUBCOVERS  = -1;
   cpxControl.IMPLED     = -1;
   cpxControl.MIRCUTS    = -1;
end
   
cpxControl

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
     cplex(c, A, x_L, x_U, b_L, b_U,  cpxControl, ...
            callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2);
toc

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / Cplex solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');

   fprintf('Inform = %d. ',Inform);
   [ExitText,ExitFlag] = cplexStatus(Inform);
   fprintf(ExitText);
   fprintf('\n');
       
   fprintf('\nObjective function at x (obj) %25.16f --- ',f_k);
   if MIP
      fprintf('Nodes visited%7d. ',glnodes);
   else
      fprintf('LP iterations%7d. ',lpiter);
   end
   fprintf('\n');
end

if PriLev > 1
   if isempty(MAX_x)
      MAX_x=length(x);
   end
   fprintf('Optimal x = \n');
   xprinte(x(1:min(length(x),MAX_x)),'x:  ');
end

if PriLev > 2
   fprintf('Slack variables s =\n');
   xprint(slack,'s:');
end

if PriLev > 3
   if isempty(MAX_c)
      MAX_c=20;
   end
   fprintf('Dual variables (Lagrangian multipliers) v = \n');
   xprinte(v(1:min(length(v),MAX_c)),'v:');

   fprintf('Reduced costs r =\n');
   xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'r: ',' %14.9f',5);
end
if PriLev > 4
   fprintf('Basis b =\n');
   xprint(basis(1:min(length(basis),MAX_x)),'b: ',' %14.9f',5);
end

% MODIFICATION LOG:
%
% 020922 hkh  Revised from Xpress to Cplex
% 020922 hkh  Adding input parameter DEFPARAM
% 050113 med  Removed commented code
% 060131 ango Revise some comments
% 070221 hkh  Revise IntVars format
% 070223 hkh  Use cplexStatus to get ExitText from Inform value
