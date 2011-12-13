% function f = qp_fX(x, Prob)
%
% Compute function value to quadratic problem P, only quadratic part
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).  f(x) = 0.5 * x'*F*x
%
% F*x is stored in global QP_Fx with corresponding x in QP_x, used by qp_g

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written June 28, 2002.     Last modified June 28, 2002.

function f = qp_fX(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic objective function

if isempty(x)
   f     = Inf;
else
   QP_x  = x;
   if isempty(Prob.QP.F)
      f = 0;
   else
      % QP Problem
      QP_Fx = Prob.QP.F*x;
      f = 0.5 * (x'*QP_Fx);
   end
end