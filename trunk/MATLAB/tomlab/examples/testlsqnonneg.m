function testlsqnonneg
echo on
% Test from OPT TB 2.0 manual page 4.110-4.111 for linear LS.
% Just use the same matrices for x>=0 problem

echo off

% Set number of calls
N = 1000;

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

fprintf('\n\n-----------------------\n\n');
fprintf('Call lsqnonneg\n\n');
fprintf('-----------------------\n\n');
fprintf('Make %d calls to see the CPU time differences\n\n',N)
tic
for i = 1:N
   [x,ResNorm,Residual,ExitFlag,Output,Lambda] = lsqnonneg(C,d,[],[]);
end
toc
fprintf('\n');

xprinte(Residual,'Residual:');
xprinte(Lambda,  'Lambda:  ');
xprinte(x,       'x:       ');
format compact
disp('Output Structure')
disp(Output)
fprintf('ResNorm %40.20f. ExitFlag %d\n',ResNorm,ExitFlag);

% Make dummy call to make the license check
Tlsqnonneg(C,d,[],[]);

fprintf('\n\n----------------------------------------------------------\n\n');
fprintf('Call Tlsqnonneg, Tomlab lsqnonneg calling mex solver Tnnls\n\n');
fprintf('----------------------------------------------------------\n\n');
Tlsqnonneg(C,d,[],[]);
fprintf('Make %d calls to see the CPU time differences\n\n',N)

tic
for i = 1:N
   [x,ResNorm,Residual,ExitFlag,Output,Lambda] = Tlsqnonneg(C,d,[],[]);
end
toc
fprintf('\n');

xprinte(Residual,'Residual:');
xprinte(Lambda,  'Lambda:  ');
xprinte(x,       'x:       ');
format compact
disp('Output Structure')
disp(Output)
fprintf('\n');