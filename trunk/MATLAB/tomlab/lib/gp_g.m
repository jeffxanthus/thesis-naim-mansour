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

function g = gp_g(x, Prob)

x=x(:);
c=Prob.GP.coef;
% Evaluate geomatric gradient

if isempty(x)
   g = zeros(Prob.N,1);
else
   g = zeros(Prob.N,1);
   f = zeros(Prob.N,1);
   for i=1:Prob.GP.not
      g(i,1) = (1/x(i)/c(i))^x(i)*log(1/c(i)/x(i))*(-1/c(i)/x(i)^2);
      g(i,1) = g(i,1);
   end
   
   f = prod((1./(x.*Prob.GP.coef)).^x);
   extra = sum(x(end-sum(Prob.GP.nct):end));
   f = f * extra^extra;
end

% MODIFICATION LOG
%
% 050509  med  Created