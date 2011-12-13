% secUpdat computes the least change formulation of secant updates.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 27, 1997. Last modified June 22, 1999.

function [A]=secUpdat(s,y,B,v)

vTs=v'*s;
z=(y-B*s)/vTs;
A=z*v' + v*z' - (z'*s/vTs)*v*v';