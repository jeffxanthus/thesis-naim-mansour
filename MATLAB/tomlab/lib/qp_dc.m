% function dc = qp_dc(x, Prob)
%
% Compute constraint Jacobian to problem P
%
% x      Point x where dc(x) is evaluated
% Prob   Problem structure
% dc     Constraint Jacobian, dc(x).  dc(x) = 2*Q*x + a;

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1996-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 22, 2004.   Last modified May 22, 2004.

function dc = qp_dc(x, Prob)

global QP_Qx QP_Qxx

x=x(:);

% Evaluate quadratic constraint Jacobian

nQC = size(Prob.QP.qc,2);

if ~isempty(QP_Qx) & length(x)==length(QP_Qxx)
    if all(QP_Qxx==x)
        for i=1:nQC
            dc(i,:) = QP_Qx(:,i)' + Prob.QP.qc(i).a';            
        end
    else
        for i=1:nQC
            dc(i,:) = (Prob.QP.qc(i).Q*x)'+Prob.QP.qc(i).a';
        end        
    end
else
    for i=1:nQC
        dc(i,:) = (Prob.QP.qc(i).Q*x)'+Prob.QP.qc(i).a';
    end
end

dc=dc*2;

% MODIFICATION LOG
%
% 040521  med  Written