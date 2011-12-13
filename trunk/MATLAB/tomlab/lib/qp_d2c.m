% function d2c = qp_d2c(x, Prob)
%
% Compute d2c to quadratic constraint
%
% x      Point x where d2c(x) is evaluated
% Prob   Problem structure
% d2c    d2c matrix, d2c(x) = lam'*Q*2

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1996-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 24, 2004.   Last modified May 24, 2004.

function d2c = qp_d2c(x, lam, Prob)

nQC = size(Prob.QP.qc,2);
d2c = sparse(Prob.N,Prob.N);
for i=1:nQC
    if lam(i) ~= 0
        d2c = d2c+lam(i)*Prob.QP.qc(i).Q*2;
    end
end

% MODIFICATION LOG
%
% 040521  med  Written