function testquadprog

% Simple test problem for QP

% Here we define the different items
% We could of course directly define them in the call to qpAssign instead.

% The matrix F in the quadratic term 0.5 x' F x
F = [2 0 0;0 2 0;0 0 2];
% The vector c in the linear term 0.5 c' x
c = zeros(3,1);

% The linear constraints b_L <= A x <= b_U
% You can define your inequalities in the most general way
A = [1 2 -1;1 -1 1];
b_L = [4 -2]';
b_U = b_L;

% Lower and upper bounds on x, x_L <= x <= x_U
x_L=[-10 -10 -10]';
x_U=[10 10 10]';
% An optional starting point
x_0=[0 0 0]';
% An optional lower bound on the objective function.
% Only used when running conSolve. Could be skipped.
f_Low=-1E5;

Name='Fletcher EQP pg 231'; % Optional, set the name of the problem

Prob = qpAssign(F,c,A,b_L,b_U,x_L,x_U,x_0,Name, [], [], f_Low);


% If we want to use graphical tools to display the function we can
% define the default plotting window in x

Prob.x_min=[0 0 -1]';
Prob.x_max=[2 2 1]';


disp(' ');
disp('testquadprog: Run TOMLAB qpSolve on simple test problem');

Result   = qpSolve(Prob);
PrintResult(Result);

disp('testquadprog: Solve QP using quadprog (optim TB 2.0 interface)'); 
disp(' '); 

disp('First test. See comments in the code')
% If A has both linear inequalities (only upper bounded) 
% and equalities we can do the following calls

ix = Prob.b_L==Prob.b_U;
E  = find(ix);
I  = find(~ix);

if xnargin('quadprog') >= 11
   [x,fVal,ExitFlag,Out,lambda] = quadprog(Prob.QP.F, Prob.QP.c,...
           Prob.A(I,:),Prob.b_U(I),...
           Prob.A(E,:),Prob.b_U(E),...
           Prob.x_L,Prob.x_U,Prob.x_0,[],Prob);
else
   [x,fVal,ExitFlag,Out,lambda] = quadprog(Prob.QP.F, Prob.QP.c,...
           Prob.A(I,:),Prob.b_U(I),...
           Prob.A(E,:),Prob.b_U(E),...
           Prob.x_L,Prob.x_U,Prob.x_0);
end
x

fprintf('\n\n')
disp('Second test. Exactly the same. Used transformation routine')
fprintf('\n')

% If inequalities are totally general, with both lower bounds not Inf
% and upper bounds, a TOMLAB utility solves this problem:

[AA,bb,meq]=cpTransf(Prob,2,0);

mA = size(AA,1);

if xnargin('quadprog') >= 11
   [x,fVal,ExitFlag,Out,lambda] = quadprog(Prob.QP.F, Prob.QP.c,...
           AA(meq+1:mA,:),bb(meq+1:mA),AA(1:meq,:),bb(1:meq),...
           Prob.x_L,Prob.x_U,Prob.x_0,[],Prob);
else
   [x,fVal,ExitFlag,Out,lambda] = quadprog(Prob.QP.F, Prob.QP.c,...
           AA(meq+1:mA,:),bb(meq+1:mA),AA(1:meq,:),bb(1:meq),...
           Prob.x_L,Prob.x_U,Prob.x_0);
end

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

fprintf('\n\n')
disp('Third test. Call with exact quadprog format ')
disp('Set solver qpopt using global structure otxProb ')
fprintf('\n')
global otxProb
otxProb = ProbDef;
otxProb.SolverQP = 'qpopt';
   [x,fVal,ExitFlag,Out,lambda] = quadprog(Prob.QP.F, Prob.QP.c,...
           AA(meq+1:mA,:),bb(meq+1:mA),AA(1:meq,:),bb(1:meq),...
           Prob.x_L,Prob.x_U,Prob.x_0);
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
otxProb = [];
