% expQG is a small example problem for defining and solving exponential
% fitting problems using the TOMLAB format. See the TOMLAB manual for
% options.

% f(t) = alpha(1)*exp(-beta(1)*t) + alpha(2)*exp(beta(2)*t)

t = [0 1.00 2.00 4.00 6.00 8.00 10.00 15.00 20.00]';
y = [905.10 620.36 270.17 154.68 106.74 80.92 69.98 62.50 56.29]';
p      = 2;                         % Two terms
Name   = 'Simple two-term exp fit'; % Problem name, can be anything
wType  = 0;                         % No weighting
SepAlg = 0;                         % Separable problem
Prob = expAssign(p,Name,t,y,wType,[],SepAlg);

Prob.SolverL2 = 'nlssol';
Result = tomRun('expSolve',Prob,1);
%Prob.SolverL2 = 'snopt';
%Result = tomRun('expSolve',Prob,1);