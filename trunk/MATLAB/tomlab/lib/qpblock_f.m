% function f = qpblock_f(x, Prob)
%
% Compute function value to quadratic factorized problem
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).
%
%           Case 1:
%
%           f(x) = 0.5 * x' * F * x + d' * x +
%           0.5 * x' * Fb.out' * Fb.inn * Fb.out * x (for i=1...p)
%
%           Case 2:
%
%           f(x) = 0.5 * x' * F * x + d' * x +
%               0.5 * x' * Fb.out' * Fb.out * x (for i=1...p)
%
% The derivative of f(x) is stored in global QP_Fx with corresponding
% x in QP_x, used by qp_g
%
%           g(x) = F * x + d  + Fb.out' * Fb.inn * Fb.out * x (for i=1...p)
%           or
%           g(x) = F * x + d  + Fb.out' * Fb.out * x (for i=1...p)

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written Aug 21, 2006.  Last modified Aug 22, 2006.

function f = qpblock_f(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic objective function

if Prob.QP.case == 1 | Prob.QP.case == 3
    if isempty(Prob.QP.F)
        % LP Problem
        if ~isempty(Prob.QP.c)
            QP_Fx = Prob.QP.c;
        else
            QP_Fx = zeros(length(x),1);
        end
    else
        % QP Problem
        if Prob.QP.case == 3
            QP_Fx = Prob.QP.F.*x;
        else
            QP_Fx = Prob.QP.F*x;
        end
        if ~isempty(Prob.QP.c)
            QP_Fx = QP_Fx + Prob.QP.c;
        end
    end
    for i=1:Prob.QP.mQP
        QP_Fx = QP_Fx + Prob.QP.Fb(i).out'*(Prob.QP.Fb(i).inn * (Prob.QP.Fb(i).out * x));
    end
else
    if isempty(Prob.QP.F)
        % LP Problem
        if ~isempty(Prob.QP.c)
            QP_Fx = Prob.QP.c;
        else
            QP_Fx = zeros(length(x),1);
        end
    else
        % QP Problem
        if Prob.QP.case == 4
            QP_Fx = Prob.QP.F.*x;
        else
            QP_Fx = Prob.QP.F*x;
        end
        if ~isempty(Prob.QP.c)
            QP_Fx = QP_Fx + Prob.QP.c;
        end
    end
    for i=1:Prob.QP.mQP
        QP_Fx = QP_Fx + Prob.QP.Fb(i).out'*(Prob.QP.Fb(i).out * x);
    end
end
if isempty(Prob.QP.c)
    f = 0.5*(x'*QP_Fx);
else
    f = 0.5*(x'*(QP_Fx+Prob.QP.c));
end
QP_x  = x;

% MODIFICATION LOG
%
% 060821  med  Written, based on qp_f
% 060822  hkh  Create derivative in QP_Fx, avoid multiplications, clean up