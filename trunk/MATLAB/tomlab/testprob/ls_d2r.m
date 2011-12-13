% function d2r = ls_d2r(x, r, J, Prob)
%
% Computes the 2nd part of the second derivative to the nonlinear 
% least squares problem at the point x for the test problem Prob.P

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function d2r = ls_d2r(x, r, J, Prob)

P=Prob.P;
uP=Prob.uP;
d2r=[];

if P==1
   % Powell
   d2r= r(2)*2*uP(1);
end

% MODIFICATION LOG:
%
% 981023  hkh  Fix 2nd derivative only for P=1
% 981108  hkh  Use general global variable LS_A instead of circle_t
% 080603  med  Switched to clsAssign, cleaned