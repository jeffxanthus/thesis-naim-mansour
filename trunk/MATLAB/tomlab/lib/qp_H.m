% function H = qp_H(x, Prob)
%
% Compute Hessian to quadratic problem
%
% x      Point x where H(x) is evaluated
% Prob   Problem structure
% H      Hessian matrix, H(x) = F, in f(x) = 0.5 * x'*F*x + c'*x;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1996-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written Sept 10, 1996.   Last modified Aug 22, 2006.

function H = qp_H(x, Prob)
H = Prob.QP.F;

% MODIFICATION LOG
%
% 990825  hkh  Modifications
% 020928  hkh  Avoid x(:), use sparse if big problem
% 060822  med  All checks moved to assign routines