% soc_dc.m
%
% function dc=soc_dc(x, Prob, varargin)
%
% soc_dc computes the constraint gradient for the soc constraints

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written Dec 22, 2000. Last modified Jul 17, 2009.

function dc=soc_dc(x, Prob, varargin)

F=Prob.SOC.A;
b=Prob.SOC.b;
C=Prob.SOC.C;
d=Prob.SOC.d;
N=Prob.SOC.N;

M = sum(N);
m = length(N);
n = length(x);
dc = zeros(m,n);

i1=1;
for i = 1:m
    i2    = i1 + N(i) -1;
    v     = F(i1:i2,:)*x - b(i1:i2);
    vNorm = norm(v);
    if vNorm == 0
       dc(i,:) = - C(i,:);
    else
       dc(i,:) = [1/norm(v)*(F(i1:i2,:)'*v)]' - C(i,:);
    end
    i1   = i2 + 1;
end

% MODIFICATION LOG
%
% 090717 med  f_0 calculation updated
