% Simple test function for fminunc, example page 4-62 in OPTB 2.0 manual
%
% function [f,g]=fminunc_fg(x)
%
% f = sin(x)+3;
% g = cos(x);

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 29, 1999.   Last modified July 29, 1999.

function [f,g]=fminunc_fg(x)

f = sin(x)+3;
g = cos(x);