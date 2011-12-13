% mpecDemo.m 
%
% Demonstrates how to setup and solve a problem with equilibrium
% (complementary) constraints. 
%
% TOMLAB /KNITRO is currently the only solver capable of solving general 
% complementarity problems 
%
% The problem is the following: 
%
% minimize  f = (x(1)-5)^2 + (2*x(2)+1)^2
%
% subject to     
%
%    2*(x(2)-1) - 1.5*x(1) + x(3) - 0.5*x(4) + x(5) == 0   % c1
%     3*x(1) - x(2) - 3     >= 0                           % c2
%      -x(1) + 0.5*x(2) + 4 >= 0                           % c3
%      -x(1) -     x(2) + 7 >= 0                           % c4
%
%    x(1:5) >= 0
%
% All four constraints are in fact linear, but are initially modelled as nonlinear
% constraints with the constants moved to the LHS.
%
% A second example will show how to code the same problem with linear constraints. 
% 
% The following complementarity demands are also imposed: 
%
% 0 <= c2(x) _|_ x(3) >= 0
%
% 0 <= c3(x) _|_ x(4) >= 0
%
% 0 <= c4(x) _|_ x(5) >= 0
%
% I.e., three of the constraints should each be complementary to a
% variable. This means that either the constraint or the variable (or both) should be fixed 
% to zero, but both can not be nonzero at the same time. 
%

% No (explicitly) linear constraints.

A   = [];
b_L = [];
b_U = [];

% Provide a pattern for the nonlinear constraints
ConsPattern = [
   1 1 1 1 1
   1 1 1 0 0
   1 1 0 0 0
   1 1 0 0 0
   ];

HessPattern = [];

% First constraint is equality == 0, the remaining are >= 0
c_L = [0,0,0,0]';
c_U = [0,inf,inf,inf]';

x_L = zeros(5,1);
x_U = inf*ones(5,1);
x_0 = x_L;

Name = 'JF Bard 1998 MPEC';

% Functions for calculating the nonlinear function and derivative values.
f   = 'mpd_f';
g   = 'mpd_g';
H   = 'mpd_H';
c   = 'mpd_c';
dc  = 'mpd_dc';
d2c = 'mpd_d2c';

Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
   [],[], ...
   A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U);

Prob.PriLevOpt = 2;

% KNITRO options. Algorithm 3 works good on this problem. 
opts     = [];
opts.ALG = 3;
Prob.KNITRO.options = opts;

% The problem structure Prob at this point describes the nonlinear program
% problem outlined earlier *without* any complementarity properties. To 
% have this automatically handled, the user now needs to specify the
% complementarity pairs:

mpec = [ ...
   3,0,0,0,0,2,0; ...
   4,0,0,0,0,3,0; ...
   5,0,0,0,0,4,0; ...
   ];

% Each row of mpec is one pair. The first says x(3) _|_ c(2) should be
% complementary. Exactly two nonzeros per row is allowed.

% The following call will return a *new* problem structure where slack
% variables have been added to handle any constraints that are part of a
% complementarity pair. 

Prob = BuildMPEC(Prob,mpec);

% Now solve this problem with KNITRO: 
Result = tomRun('knitro',Prob,2);

% The original problem is available as Prob.orgProb.

N0 = Prob.orgProb.N;

x = Result.x_k;
s = x(N0+1:end)     % The slacks added by BuildMPEC
x = x(1:N0)         % The "original" variables

% Constraint values for the modified problem. If not infeasible, these
% should all be zero because BuildMPEC has changed the constraints in the
% mpec pairs to include slack variables. 

c = Result.c_k           

% The permutation matrix for the slacks in the nonlinear constraints is
% available: 

MP = full(Prob.MPEC.MP);

% The values of the constraints + their respective slacks:

C = c + MP*s

x(3:5), C(2:4) 

% This should give something like: 
%
% ans =
% 
%     3.5000
%          0
%          0
% 
% 
% ans =
% 
%      0
%      3
%      6
%
% I.e, we have the complementarity conditions satisfied. 

pause

% Save the results we have so far. 

Prob1   = Prob;
Result1 = Result;

clear Prob Result

% Now lets look at how to do this with a linear constraints. 
% In fact, this is a quadratic problem (QP), so let's set it up as one: 

% minimize  f = (x(1)-5)^2 + (2*x(2)+1)^2
%
% subject to     
%
%    2*(x(2)-1) - 1.5*x(1) + x(3) - 0.5*x(4) + x(5) == 0   % c1
%     3*x(1) - x(2) - 3     >= 0                           % c2
%      -x(1) + 0.5*x(2) + 4 >= 0                           % c3
%      -x(1) -     x(2) + 7 >= 0                           % c4
%
%    x(1:5) >= 0

% We can write f = x'*F*x + c'*x (apart from a missing constant term 27:

F = zeros(5,5);
F(1,1) = 1;
F(2,2) = 4;
c      = [-10,4,0,0,0]';

% Linear constraints, by moving constant terms in c1-c4 to the right hand
% side: 
A = [...
   -1.5   2.0  1.0  -0.5  1.0 ; ...
    3.0  -1.0  0     0    0   ; ...
   -1.0   0.5  0     0    0   ; ...
   -1.0  -1.0  0     0    0 ];

b_L = [ 2 , 3   , -4  , -7 ]';
b_U = [ 2 , inf , inf , inf]';

% Lower and upper bounds
x_L =    zeros(5,1);
x_U = inf*ones(5,1);
x_0 = [];

Name = 'JF-BARD-1998-QP';

% Assign a TOMLAB 'qp' problem: 
Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name);

disp('We still want c2 _|_ x3, c3 _|_ x4, c4 _|_ x5, but now the constraints are linear: ')

mpec = [ ...
      3,0, 2,0, 0,0 ; ...
      4,0, 3,0, 0,0 ; ...
      5,0, 4,0, 0,0 ]

Prob = BuildMPEC(Prob,mpec);

disp('Three slacks have been added to the problem, easy to see by looking at the new linear constraint matrix: ');
A = Prob.A

disp('... and the old one:')
A_orig = Prob.orgProb.A

pause

% Solve the QP (with MPEC pairs) using KNITRO:
Prob.PriLevOpt = 2;

% Enable some crossover iterations, to "polish" the solution a bit in case KNITRO chooses an interior point algorithm: 
Prob.KNITRO.options.MAXCROSSIT = 100;

Result = tomRun('knitro',Prob,2);

x = Result.x_k

% Linear constraints result:
Ax = Result.Ax

% The slacks..
s = x(6:8)
