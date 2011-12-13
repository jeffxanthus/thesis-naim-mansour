% function f = ego_c(x, Prob)
%
% ego_c does rescaling of x, if Prob.SCALE is true, 
% before calling the user constraint routine.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2008 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Apr 4, 2002.   Last modified Jan 29, 2008.

function c = ego_c(x, Prob)

if Prob.cNargin > 0
   if Prob.SCALE > 0
      x = tomsol(9, Prob.xL, x ,Prob.xD); 
   end

   if Prob.cNargin==2
      c=feval(Prob.c, x, Prob);
   else
      c=feval(Prob.c, x);
   end
   c = c(:);
else
   c = [];
end

% MODIFICATION LOG:
%
% 020404  hkh  Written
% 021020  hkh  Avoid case when Prob.c undefined, check on Prob.cNargin
% 080129  hkh  Avoid additional constraint, only user constraints
