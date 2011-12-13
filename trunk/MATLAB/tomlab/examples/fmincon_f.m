% Simple test function for fmincon, example page 4-37 in OPTB 2.0 manual
%
% function f=fmincon_f(x)
%
% f= -x(1) * x(2) * x(3);

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 29, 1999.   Last modified July 29, 1999.

function f=fmincon_f(x)

f= -x(1) * x(2) * x(3);