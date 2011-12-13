% cpxiptest.m:
%
% Test of TOMLAB /CPLEX MEX interface solving three integer linear 
% optimization problems calling CPLEX solver
%
% The test problems have 61 variables and 138 linear inequalities
% 32 of the 138 inequalities are just zero rows in the matrix A.
%
% The three problems are stored in ilp061.mat, ilp062.mat and ilp063.mat.
%
% function cpxiptest(Cut, PreSolve, cpxControl)
%
% The input determines the change of five different control parameters 
%
% Cut       Value of the cut strategy control parameters, Default Cut = -1
%           -1 = Auto select of 1 or 2, 0 = No cuts, 1 = Conservative cuts, 
%           2  = Aggressive cut strategy
%
% PreSolve  Value of the PREIND control parameter
%           0 = No presolve, 1 = presolve. Default = 1
%
% cpxControl The initial cpxControl structure. Here the user may set 
%            additional control parameters

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Sept 27, 2002. Last modified Feb 23, 2007.

function cpxiptest(Cut, PreSolve, cpxControl)

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


fprintf('\n\n*******************************************\n');
fprintf('cpxiptest: set up problem 1:\n');

load('ilp061.mat');
if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end


if 1
   Prob.P = 1;
   x_U = NewxU(x_U,A,b_U);
   xprinti(x_U,'x_U');

   % Find 0 rows and remove them
   ix = find(sum(A' > 0) ~= 0);
   % xprinti(ix,'ix0');
   ix0 = find(sum(A' > 0) == 0);
   rows_0=length(ix0)

   b_L = b_L(ix);
   b_U = b_U(ix);
   A   = A(ix,:);
end


disp('cpxiptest: Solving problem 1...')

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

xprinti(x,'x:')
maxx=max(x)

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);


fprintf('\n\n*******************************************\n');
fprintf('cpxiptest: set up problem 2:\n');

load('ilp062.mat');
if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end


if 1
   Prob.P = 2;
   x_U = NewxU(x_U,A,b_U);
   xprinti(x_U,'x_U');

   % Find 0 rows and remove them
   ix = find(sum(A' > 0) ~= 0);
   % xprinti(ix,'ix0');
   ix0 = find(sum(A' > 0) == 0);
   rows_0=length(ix0)
   b_L = b_L(ix);
   b_U = b_U(ix);
   A   = A(ix,:);
end

disp('cpxiptest: Solving problem 2...')

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

xprinti(x,'x:')
maxx=max(x)
cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);

fprintf('\n\n*******************************************\n');
fprintf('cpxiptest: set up problem 3:\n');

load('ilp063.mat');
if length(noivars) == 1
   IntVars = [1:noivars];
else
   IntVars = noivars;
end


if 1
   Prob.P = 3;
   x_U = NewxU(x_U,A,b_U);
   xprinti(x_U,'x_U');

   % Find 0 rows and remove them
   ix = find(sum(A' > 0) ~= 0);
   % xprinti(ix,'ix0');
   ix0 = find(sum(A' > 0) == 0);
   rows_0=length(ix0)
   b_L = b_L(ix);
   b_U = b_U(ix);
   A   = A(ix,:);
end

disp('cpxiptest: Solving problem 3...')

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, [], [], IntVars);
toc

xprinti(x,'x:')
maxx=max(x)

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc);

function x_U = NewxU(x_U,A,b_U)

% Squeeze the upper bounds on x
% x >= 0, all elements in A >= 0, All b_U >= 0

n = size(A,2);

for i = 1:n
    z = A(:,i);
    ix = find(z > 0);
    x_U(i) = floor(min(b_U(ix)./z(ix)));
end

% MODIFICATION LOG:
%
% 020927 hkh  Revision for CPLEX
% 021001 hkh  Add cuts and presolve
% 021003 ango Add explanatory texts
% 050209 med  Changed name from iptest to cpxiptest
% 061130 ango Trivial errors in printouts fixed
% 070223 hkh  Generate IntVars from noivars
