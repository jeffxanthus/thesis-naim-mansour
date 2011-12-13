% Test of TOMLAB /CPLEX.
%
% Running a generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
%
% Define the linear sos1 constraints explicitly
%
% function x = cpxTest1(DEFPARAM)
%
% DEFPARAM  If =1 use default parameters, presolve, cuts, Dual simplex, 
%           If =0 no presolve or cuts, Primal simplex (default)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Oct 28, 1999.   Last modified Feb 23, 2007.

function x = cpxTest1(DEFPARAM)
if nargin < 1
   DEFPARAM = [];
end

if isempty(DEFPARAM), DEFPARAM = 0; end
  
% Global Tomlab variables that determine number of elements printed
global MAX_x MAX_c


A=[ 3 24 53 27 17; ...
   15 23 43 74 23; ...
   54 43 27 21 36; ...
   92 83 45 35 23; ...
   19 10 33 43 12; ...
   91 55 32 26 23; ...
   15 25 35 37 28; ...
   47 43 35 32 37; ...
   34 23 52 46 43; ...
   35 23 34 25 40];

C=-[15 44 76 43 34; ...
   19 23 45 46 34; ...
   10  6  3 23 15; ...
   60 45 34 36 23; ...
   21 12 34 44 10; ...
   67 35 34 20 37; ...
   23 34 44 47 32; ...
   23 25 32 15 27; ...
   15 13 23 24 34; ...
   10 15 23 12 13];

b = [80 63 75 98 59]';

MIP=1;

[c,x_L,x_U,b_L,b_U,A]=abc2gap(A,b,C,0);


IntVars=[1:length(c)];
MAX_x = length(c);

% PriLev > 0 gives printing in this file
% Increasing PriLev gives more output

PriLev = 1;

callback=[]; % Use default callbacks 
PI      = [];
SC      = [];
SI      = [];
sos1    = [];
sos2    = [];

Prob.P = 1;  % Initialize Prob with problem number 1

% IntVars=[]
% MIP = 0

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
% 020922 hkh Revised from Xpress to Cplex
% 020922 hkh Adding input parameter DEFPARAM
% 050113 med Removed commented code
% 050117 med mlint review
% 070221 hkh Revise IntVars format
% 070223 hkh Get ExitText from Inform using cplexStatus
