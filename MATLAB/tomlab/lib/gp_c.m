% function c = gp_c(x, Prob)
%
% Compute constraint values to geometric programming problem (primal).
%
% x      Point x where c(x) is evaluated
% Prob   Problem structure
% c      Constraint vector, c(x).

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.8.0$
% Written May 9, 2005.  Last modified May 9, 2005.

function c = gp_c(x, Prob)

x=x(:);
clen = length(Prob.GP.nterm)-1;
c = zeros(clen,1);
for i=1:clen
   nterm = Prob.GP.nterm(i+1);
   for j=1:nterm
      c(i,1) = c(i,1) + Prob.GP.coef(sum(Prob.GP.nterm(1:i))+j)*prod(x.^(Prob.GP.A(sum(Prob.GP.nterm(1:i))+j,:)'));
   end
end

% MODIFICATION LOG
%
% 060608  med  Constraints implemented on primal form