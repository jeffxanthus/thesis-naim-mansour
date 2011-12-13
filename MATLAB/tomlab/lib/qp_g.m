% function g = qp_g(x, Prob)
%
% Compute gradient to quadratic problem P
%
% x      Point x where g(x) is evaluated 
% Prob   Problem structure
% g      Gradient, g(x).  g(x) = F*x + c;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1996-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Sept 10, 1996.   Last modified Sep 5, 1999.

function g = qp_g(x, Prob)

global QP_x QP_Fx

x=x(:);

% Evaluate quadratic gradient

if isempty(x)
    g = [];
    return
end

if isequalwithequalnans(QP_x,x)
   g = QP_Fx;
else
   f = qp_f(x, Prob);
   g = QP_Fx;
end