% function g = qp_gX(x, Prob)
%
% Compute gradient to quadratic problem P, only quadratic part
%
% x      Point x where g(x) is evaluated
% Prob   Problem structure
% g      Gradient, g(x).  g(x) = F*x

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written June 28, 2002.     Last modified June 28, 2002.

function g = qp_gX(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic gradient

if isempty(x)
   g = [];
   return
end
if isempty(QP_Fx)
   g = Prob.QP.F*x;
else
   g = QP_Fx;
   %if all(QP_x==x)
   %   g = QP_Fx + c;
   %else
   %   g = Prob.QP.F*x;
   %end
end