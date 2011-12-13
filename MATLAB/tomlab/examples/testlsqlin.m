function testlsqlin
echo on
% Test from OPT TB 2.0 manual page 4.110-4.111
echo off
global otxProb
otxProb = [];

C=[ 0.9501   0.7620   0.6153   0.4057
    0.2311   0.4564   0.7919   0.9354
    0.6068   0.0185   0.9218   0.9169
    0.4859   0.8214   0.7382   0.4102
    0.8912   0.4447   0.1762   0.8936];

d = [ 0.0578 0.3528 0.8131 0.0098 0.1388]';

A =[ 0.2027   0.2721   0.7467   0.4659
     0.1987   0.1988   0.4450   0.4186
     0.6037   0.0152   0.9318   0.8462];

b =[ 0.5251 0.2026 0.6721]';

lb = -0.1*ones(4,1);
ub = 2*ones(4,1);
 
% Use if want to change solver in lsqlin
% Prob=ProbDef;
% Prob.SolverQP='lssol';
% [x,ResNorm,Residual,ExitFlag,Output,Lambda] = lsqlin(C,d,A,b,[],[],lb,ub,...
% [],[],Prob);

[x,ResNorm,Residual,ExitFlag,Output,Lambda] = lsqlin(C,d,A,b,[],[],lb,ub);

xprinte(Residual,'Residual:      ');
xprinte(A*x-b,'Constraint Ax-b');
xprinte(Lambda.eqlin,'Lambda.eqlin:  ');
xprinte(Lambda.ineqlin,'Lambda.ineqlin:');
xprinte(Lambda.lower,'Lambda.lower:  ');
xprinte(Lambda.upper,'Lambda.upper:  ');
xprinte(x,'x:             ');
format compact
disp('Output Structure')
disp(Output)
fprintf('ResNorm %40.20f. ExitFlag %d\n',ResNorm,ExitFlag);


fprintf('\n');
disp('The 2nd constraint is fulfilled.');
disp('Change the call, setting 2nd constraint as equality.');
disp('Just to test that the solution is the same.');
fprintf('\n');

if exist('optimset')
   options=optimset;
else
   options=[];
end
%options.Diagnostics='on';
%options.Display='iter';


[x,ResNorm,Residual,ExitFlag,Output,Lambda] = lsqlin(C, d, ...
   A([1 3],:),b([1 3]),A(2,:),b(2),lb,ub,[],options);

xprinte(Residual,'Residual:      ');
xprinte(A*x-b,'Constraint Ax-b');
xprinte(Lambda.eqlin,'Lambda.eqlin:  ');
xprinte(Lambda.ineqlin,'Lambda.ineqlin:');
xprinte(Lambda.lower,'Lambda.lower:  ');
xprinte(Lambda.upper,'Lambda.upper:  ');
xprinte(x,'x:             ');
format compact
disp('Output Structure')
disp(Output)
fprintf('ResNorm %40.20f. ExitFlag %d\n',ResNorm,ExitFlag);