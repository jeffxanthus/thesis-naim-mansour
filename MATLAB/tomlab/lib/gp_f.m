% function f = gp_f(x, Prob)
%
% Compute function value to geometric programming problem (primal).
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.8.0$
% Written May 9, 2005.  Last modified May 9, 2005.

function f = gp_f(x, Prob)

x=x(:);

nterm = Prob.GP.nterm(1);
f = 0;

for i=1:nterm
   f = f + Prob.GP.coef(i)*prod(x.^(Prob.GP.A(i,:)'));
end

% MODIFICATION LOG
%
% 050509  med  Created
% 050510  med  Updated help
% 050519  frhe Now returns the non-negated function value of the dual
%              problem. This function should be maximized that is.
% 060608  med  Objective implemented on primal form
