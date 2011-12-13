% function f = ego_cc(x, Prob)
%
% ego_cc does rescaling of x, if Prob.SCALE is true, 
% before calling the user constraint routine.
% It first calls ego_c with the original x

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Apr 4, 2002.   Last modified Oct 20, 2002.

function c = ego_cc(x, Prob)

c1 = ego_c(x,Prob);

if Prob.cNargin > 0
   if Prob.SCALE > 0
      x = tomsol(9, Prob.xL, x ,Prob.xD); 
   end

   if Prob.cNargin==2
      c=feval(Prob.c, x, Prob);
   else
      c=feval(Prob.c, x);
   end
   c = [c(:);c1];
else
   c = c1;
end

% MODIFICATION LOG:
%
% 020404  hkh  Written
% 021020  hkh  Avoid case when Prob.c undefined, check on Prob.cNargin