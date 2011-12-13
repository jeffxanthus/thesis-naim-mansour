% function f = qp_f(x, Prob)
%
% Compute function value to quadratic problem P
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).  f(x) = 0.5 * x'*F*x + c'*x;
%
% F*x is stored in global QP_Fx with corresponding x in QP_x, used by qp_g

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written April 28, 1995.  Last modified Aug 22, 2006.

function f = qp_f(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic objective function
QP_x  = x;
if isempty(Prob.QP.c)
    QP_Fx = Prob.QP.F*x;
    f = 0.5*(x'*QP_Fx);
else
    QP_Fx = Prob.QP.F*x+Prob.QP.c;
    f = 0.5*(x'*(QP_Fx+Prob.QP.c));
end