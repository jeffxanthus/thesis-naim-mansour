% soc_c.m
%
% function c = soc_c(x, Prob)
%
% soc_c computes the SOC constraints in the point x

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Dec 22, 2000. Last modified Dec 22, 2000.

function c = soc_c(x, Prob)

F=Prob.SOC.A;
b=Prob.SOC.b;
C=Prob.SOC.C;
d=Prob.SOC.d;
N=Prob.SOC.N;

M = sum(N);
m = length(N);
c = zeros(m,1);

i1=1;
for i = 1:m
    i2   = i1 + N(i) -1;
    v    = F(i1:i2,:)*x - b(i1:i2);
    c(i) = norm(v) - C(i,:)*x;
    i1   = i2 + 1;
end

% MODIFICATION LOG
%