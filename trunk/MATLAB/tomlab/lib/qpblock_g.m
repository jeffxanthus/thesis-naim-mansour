% function g = qpblock_g(x, Prob)
%
% Compute gradient for quadratic factorized problem
%
% x      Point x where g(x) is evaluated
% Prob   Problem structure
% g      Gradient, g(x). See qpblock_f for more information.
%
% Use global QP_Fx value from qpblock_f if QP_x == x, otherwise call
% qpblock_f to get correct value

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written Aug 21, 2006.   Last modified Aug 22, 2006.

function g = qpblock_g(x, Prob)

global QP_x QP_Fx

x=x(:);

if all(QP_x==x)
   g = QP_Fx;
else
   f = qpblock_f(x, Prob);
   g = QP_Fx;
end

% MODIFICATION LOG
%
% 060821  med  Written, based on qp_g
% 060821  hkh  Revised code, only use QP_Fx from qpblock_f, safe guarded