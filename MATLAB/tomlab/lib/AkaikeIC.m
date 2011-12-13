% function f = AkaikeIC(ssq, N, p)
%
% Compute AIC; Akaike's information criterion (approximated)
%
% ssq   Sum of squares of the residuals
% N     Number of data points
% p     Number of model parameters
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Aug 10, 2002.   Last modified Feb 11, 2005

function f = AkaikeIC(ssq, N, p)

f = N*log(ssq/N) + 2*p;