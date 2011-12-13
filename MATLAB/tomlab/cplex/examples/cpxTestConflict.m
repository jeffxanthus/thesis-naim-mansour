% Demonstration of the TOMLAB /CPLEX Conflict Refinement feature
%
% Running a generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
%
% Define the linear sos1 constraints explicitly. Modify a bound to produce
% an infeasibility and invoke CPLEX again with Conflict Refinement enabled.
%
% function x = cpxTestConflict()
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Jan 31, 2006.   Last modified Feb 21, 2007.

function x = cpxTestConflict

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

% Increasing PriLev gives more output
PriLev = 2;

callback= [];
PI      = [];
SC      = [];
SI      = [];
sos1    = [];
sos2    = [];

Prob    = [];

% CPLEX parameters structure
cpxControl = [];

F = []; % No quadratic term in objective function
logfile = '';
savefile= '';
savemode = [];

% No quadratic constraints
qc       = [];

% No Conlict File requested. 
conflictFile = '';

% Solve problem, should not pose any difficulties
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
   cplex(c, A, x_L, x_U, b_L, b_U,  cpxControl, ...
   callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2);

% Change the first constraint bounds to something that can not be satisfied: 
b_L(1)=1;
b_U(1)=1;

% Build a conflict group descriptor: 
confgrps = cpxBuildConflict( length(x_L), length(b_L), 0, 0, 0, 'full');

[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, ...
   confstat, iconfstat, sa] = ...
   cplex(c, A, x_L, x_U, b_L, b_U,  cpxControl, ...
   callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2, F, ...
   logfile, savefile, savemode, qc, confgrps, conflictFile);

% Find conflict members. These will have a confstat(k).istat value greater
% than -1 (excluded)
%
% The iconfstat array could also have been investigated to find elements that
% are greater than zero and thus potential conflict members. 
%
for k=1:length(confstat)
   if(confstat(k).istat>-1)
      c = confstat(k);
      fprintf('Conflict group %d: %s (%d)\n',k,c.status,c.istat)
      c
   end
end

% The resulting output will show that conflict group 101, consisting of the
% first linear constraint is be a member of the conflict (in fact, it is the 
% only member in the conflict).

% MODIFICATION LOG:
%
% 060131 ango Wrote file
% 060203 ango Format changes
% 060210 ango Minor printout fix
% 070221 hkh  Revise IntVars format
