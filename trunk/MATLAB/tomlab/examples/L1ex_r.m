% function r=L1ex_r(x,Prob)
%
% Compute the residual for the L1 problem Madsen-Tinglett I

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Mar 28, 2002.   Last modified Dec 14, 2002.

function r=L1ex_r(x,Prob)

r(1)=-x(1)*(1-x(2)); 
r(2)=-x(1)*(1-x(2)^2);
r(3)=-x(1)*(1-x(2)^3);

r=r(:);

r = r - Prob.LS.y;

%xprint(x,'x:');
%xprinte(r,'r:');