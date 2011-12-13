% function J=L1ex_J(x,Prob)
%
% Compute the Jacobian for the L1 problem Madsen-Tinglett I

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Mar 28, 2002.   Last modified Dec 14, 2002.

function J=L1ex_J(x,Prob)

J = [x(2)-1,x(1);x(2)^2-1,2*x(1)*x(2);x(2)^3-1,3*x(1)*x(2)^2];