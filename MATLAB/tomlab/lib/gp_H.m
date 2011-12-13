% function f = gp_f(x, Prob)
%
% Compute function value to geometric programming problem.
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).  f(x) = -f(x);

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.8.0$
% Written May 9, 2005.  Last modified May 9, 2005.

function f = gp_H(x, Prob)

x=x(:);

% Evaluate geometric gradient

if isempty(x)
   f = Inf;
else
   f = prod((1./(x.*Prob.GP.coef)).^x);
   extra = sum(x(end-sum(Prob.GP.nct):end));
   f = f * extra^extra;
end

% MODIFICATION LOG
%
% 050509  med  Created