% function c = qp_c(x, Prob)
%
% Compute quadratic constraints to problem P
%
% x      Point x where c(x) is evaluated
% Prob   Problem structure
% c      Constraint values, c(x).  c(x) = x'*Q*x + a'*x;
%
% Q*x is stored in global QP_Qx with corresponding x in QP_Qxx, used by
% qp_dc

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 21, 2004.  Last modified May 21, 2004.

function c = qp_c(x, Prob)

global QP_Qx QP_Qxx
x=x(:);
% Evaluate quadratic constraints

nQC = size(Prob.QP.qc,2);
c = zeros(nQC,1);
for i=1:nQC
    QP_Qx(:,i) = Prob.QP.qc(i).Q*x;
    c(i,1) = x'*QP_Qx(:,i) + Prob.QP.qc(i).a'*x;
end

QP_Qxx = x;

% MODIFICATION LOG
%
% 040521  med  Written