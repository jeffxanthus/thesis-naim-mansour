% Test of TOMLAB /CPLEX. Knapsack problems
%
% function cpxKnaps(P, Cut)
%
% INPUT:
% P        Problem number 1-3. Default 1.
%
% Cut      Cut strategy. 0 = no cuts, 1 = cuts, 2 = aggressive cuts
%          Default 0.
%
% Problem 1: 'Weingartner 1 - 2/28 0-1 knapsack';
% Problem 2: 'Hansen, Plateau 1 - 4/28 0-1 knapsack';
% Problem 3: 'PB 4 - 2/29 0-1 knapsack';

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written Dec 19, 1999.      Last modified Feb 23, 2007.

function cpxKnaps(P, Cut)

if nargin < 2
   Cut=[];
   if nargin < 1
      P=[];
   end
end
if isempty(P),   P=1; end
if isempty(Cut), Cut=0; end

if P < 1 | P > 3
   P=1;
end
if Cut < 0 | Cut > 2
   Cut=0;
end
  
global MAX_x MAX_c

if P==1
   Name='Weingartner 1 - 2/28 0-1 knapsack';
   A = [ 45      0     85     150     65     95     30      0    170  0 ...
         40     25     20       0      0     25      0      0     25  0 ...
         165      0     85       0      0      0      0    100  ; ...
         30     20    125       5     80     25     35     73     12  15 ...
         15     40      5      10     10     12     10      9      0  20 ...
         60     40     50      36     49     40     19    150]; 
   b_U = [600;600];  % // 2 knapsack capacities
   c   = [1898  440  22507  270  14148   3100   4650  30800   615   4975 ...
     1160   4225    510   11880    479    440    490    330    110    560 ...
     24355   2885  11748    4550    750   3720   1950  10500]';  % 28 weights
   % True optimum
   f_opt=-141278;
elseif P==2
   Name='Hansen, Plateau 1 - 4/28 0-1 knapsack';

   A = [ 40	 91	3    12     3	 18    25     1     1	  8 ...
         1	  1    49     8    21	  6	1     5     8	  1 ...
         0	 42	6     4     8	  0    10     1 ;     ...
        16	 92	4    18     6	  0	8     2     1	  6 ...
         2	  1    70     9    22	  4	1     5     6	  0 ...
         4	  8	4     3     0	 10	0     6 ;     ...
        38	 39	5    40     8	 12    15     0     1	 20 ...
         3	  0    40     6     8	  0	6     4     4	  1 ...
         5	  8	2     8     0	 20	0     0 ;     ...
        38	 52	7    20     0	  3	4     1     2	  4 ...
         6	  1    18    15    38	 10	4     8     0	  3 ...
         0	  6	1     3     0	  3	5     4 ];
   b_U = [219	203   208   180]';
   c   = [ 560  1125    68   328    47	122   196    41    25	115 ...
           82	 22   631   132   420	 86    42   103    81	 26 ...
           49	316    72    71    49	108   116    90]';
   % True optimum
   f_opt=-3418;
elseif P==3
   Name='PB 4 - 2/29 0-1 knapsack';
   A = [ 25  17  20  22 15 10  50 10 150  0  0 0 40 60 , zeros(1,15); ...
         zeros(1,5), 2 5 6 40 2	6 10 13 30 15 5 5 10 15	91 24 15 ...
         15 5 10 15 10 10 10];			   
   b_U = [ 153	 154]';
   c   = [ 7074	5587  5500   3450   367  4235	9262  6155  32240   1600 ...
           4500	6570  7010  16020  8000  2568	2365  4350   4975  29400 ...
           7471	3840  3575   1125  1790  2450	 540   445    220]';
   % True optimum
   f_opt=-95618;

end

disp(Name)

% Change sign to make minimum problem
c=-c;

[m,n]=size(A);

% Lower and upper bounds
x_L = zeros(n,1);
x_U = ones(n,1);
b_L = -Inf*ones(m,1);


IntVars = [1:length(c)];
MIP     = 1;

MAX_x = length(c);

PriLev = 1

callback=[];

PI      = [];
SC      = [];
SI      = [];
sos1    = [];
sos2    = [];

Prob.P=P;

if Cut == 0, Cut = -1; end
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

cpxControl

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
     cplex(c, A, x_L, x_U, b_L, b_U,  cpxControl, ...
            callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2);
toc

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / Cplex solving Problem %d.',Prob.P);
   fprintf(' - %s',Name);
   fprintf('\n');
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');

   [ExitText,ExitFlag] = cplexStatus(Inform);

   fprintf('Inform = %d. ',Inform);
   fprintf(ExitText);
   fprintf('\n');
   
   fprintf('\nObjective function at x (obj) %25.16f\n',f_k);
   if MIP
      fprintf('\nTrue optimum       at x (obj) %25.16f\n\n',f_opt);
      fprintf('Nodes visited%7d. ',glnodes);
   else
      fprintf('\nLP iterations%7d. ',lpiter);
   end
   fprintf('\n');
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
% 020107 hkh Last revised for Xpress
% 020831 fre Made it work for CPLEX
% 020921 hkh Removed callback heuristic
% 050113 med Removed commented code
% 070221 hkh Revise IntVars format
% 070223 hkh Use new cplexStatus to get ExitText from cplexStatus
