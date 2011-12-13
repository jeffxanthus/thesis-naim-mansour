% function [H,Hx] = qp_Hess(x, Prob)
%
% Compute Hessian and Hessian times vector, for quadratic problem
%
% x      Point x where H(x) is evaluated
% Prob   Problem structure
% H      Hessian matrix, H(x) = F, in f(x) = 0.5 * x'*F*x + c'*x;
% Hx     H*x

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Oct 4, 1998.     Last modified Sep 16, 1999.

function [H,Hx] = qp_Hess(x, Prob)

global QP_x QP_Fx

x=x(:);

if isempty(Prob.QP.F)
   H  = zeros(length(x),length(x));
   Hx = zeros(length(x),1);
else
   H  = Prob.QP.F;
   if all(x==QP_x)
      Hx=QP_Fx;
   else
      Hx = H*x;
   end
end