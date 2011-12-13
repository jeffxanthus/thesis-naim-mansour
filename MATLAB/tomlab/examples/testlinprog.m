function testlinprog
% Simple LP problem with two variables on inequality form 
%
% min c'x subject to Ax <= b, x >= 0
%
% Here we define the different items
% We could of course directly define them in the call to lpAssign instead.

disp('Test inequality form of simple LP')

% The vector c in the linear term 0.5 c' x in the objective function

c = [-7;-5];

% The linear constraints b_L <= A x <= b_U
% You can define your inequalities in the most general way
% See testquadprog.m and qpAssign.m for more general discussion
A = [1 2 ;4 1 ];
b_U = [6 12]';
% Setting lower bounds to empty equivalent to setting them to b_L=[-Inf;-Inf];
b_L = [-Inf;-Inf];
% b_L not equal to b_U implies inequalities

% Lower and upper bounds on x, x_L <= x <= x_U
x_L=zeros(2,1);
% Upper bounds empty. Equivalent to setting them to Inf (x_U=[Inf;Inf;Inf;Inf])
x_U=[];

% An optional starting point. Set empty now
x_0=[0;0];

% An optional lower bound on the objective function.
% Only used when running conSolve. Could be skipped.
f_Low=-1E5;

% Optional, set the name of the problem
Name='Simple LP problem on inequality form'; 
Prob = lpAssign(c,A,b_L,b_U,x_L,x_U,x_0,Name,[],[],f_Low);


% If we want to use graphical tools to display the function we can
% define the default plotting window in x

Prob.x_min=[0 0]';
Prob.x_max=[4 4]';

disp(' ');
disp('testlinprog: Run TOMLAB lpSimplex on simple inequality LP test problem');

Result   = lpSimplex(Prob);
PrintResult(Result);

disp('testlinprog: Run linprog on simple inequality LP ');


[x,fVal,ExitFlag,Out,lambda]=linprog(Prob.QP.c,Prob.A,Prob.b_U,[],[],Prob.x_L);
format long
x
fVal
ExitFlag
Out
disp('lambda.ineqlin')
xprinte(lambda.ineqlin)
disp('lambda.eqlin')
xprinte(lambda.eqlin)
disp('lambda.lower')
xprinte(lambda.lower)
disp('lambda.upper')
xprinte(lambda.upper)

% Simple LP problem on Standard form (as equalities)
%
% min c'x subject to Ax=b, x >= 0

% Here we define the different items
% We could of course directly define them in the call to lpAssign instead.

% The vector c in the linear term 0.5 c' x in the objective function
% Two slack variables are added, zero coefficient in the objective
c = [-7; -5; 0; 0];

% The linear constraints b_L <= A x <= b_U. 
% A is expanded with an indentity matrix for the slack variables
A = [1 2 1 0;4 1 0 1 ];
b_L = [6 12]';
% Setting lower bounds equal to make it equalities
b_U = b_L;

% Lower and upper bounds on x, x_L <= x <= x_U
x_L=zeros(4,1);
% Upper bounds empty. Equivalent to setting them to Inf (x_U=[Inf;Inf;Inf;Inf])
x_U=[];

% An optional starting point. Set empty now
x_0=[];

% An optional lower bound on the objective function.
% Only used when running conSolve. Could be skipped.
f_Low=-1E5;

Name='Simple LP problem'; % Optional, set the name of the problem
Prob = lpAssign(c,A,b_L,b_U,x_L,x_U,x_0,Name,[],[],f_Low);


% If we want to use graphical tools to display the function we can
% define the default plotting window in x

Prob.x_min=[0 0 0 0]';
Prob.x_max=[4 4 4 4]';

disp(' ');
fprintf('testlinprog: Run TOMLAB lpSimplex on simple LP test problem ');
fprintf('on standard form, as equalities \n');

Result   = lpSimplex(Prob);
PrintResult(Result);

disp('testlinprog: Run linprog on LP on standard form, as equalities');

[x,fVal,ExitFlag,Out,lambda]=linprog(Prob.QP.c,[],[],Prob.A,Prob.b_U,Prob.x_L);
format long
x
fVal
ExitFlag
Out
disp('lambda.ineqlin')
xprinte(lambda.ineqlin)
disp('lambda.eqlin')
xprinte(lambda.eqlin)
disp('lambda.lower')
xprinte(lambda.lower)
disp('lambda.upper')
xprinte(lambda.upper)