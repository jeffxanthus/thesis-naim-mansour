% Simple test of QP problem running CPLEX
%
% function cpxTestQP1(MIP, DEFPARAM)
%
% if MIP == 1, try to make x(3) integer
% if MIP == 0, solve as a standard QP
%
% DEFPARAM  If =1 use default parameters, presolve, cuts, Dual simplex, 
%           If =0 no presolve or cuts, Primal simplex (default)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written July 10, 2001.  Last modified Feb 23, 2007.

function cpxTestQP1(MIP, DEFPARAM)

if nargin < 2
   DEFPARAM = [];
   if nargin < 1
      MIP = [];
   end
end

if isempty(DEFPARAM), DEFPARAM = 0; end
if isempty(MIP),      MIP = 0; end

% Global Tomlab variables that determine number of elements printed
global MAX_x MAX_c
MAX_x=3; MAX_c=2;

c    = zeros(3,1);
Name = 'Fletcher EQP pg 231';
F    = [2 0 0;0 2 0;0 0 2];
A    = [1 2 -1;1 -1 1];
b_L  = [4 -2]';
b_U  = b_L;

x_L  = [-10 -10 -10]';
x_U  = [10 10 10]';

callback  = [];
PriLev    = 1;
Prob.P    = 1;
PI        = [];
SC        = [];
SI        = [];
sos1      = [];
sos2      = [];

% CPLEX will return the QP solution.
if MIP
   IntVars   = logical([0 0 1]);
else
   IntVars   = [];
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
            callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2,F);
toc

PriLev = 1;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / CPLEX solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');

   fprintf('Inform = %d. ',Inform);
   [ExitText,ExitFlag] = cplexStatus(Inform);
   fprintf(ExitText);
   fprintf('\n');   
   fprintf('\nObjective function at x (obj) %25.16f --- ',f_k);
   if MIP
      fprintf('Nodes visited%7d. ',glnodes);
   else
      fprintf('QP iterations%7d. ',lpiter);
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
% 020922 hkh Revised from Xpress to Cplex
% 020922 hkh Adding input parameter DEFPARAM
% 050113 med Removed commented code
% 070221 hkh Revise IntVars format
% 070223 hkh Use cplexStatus to get ExitText from Inform value
