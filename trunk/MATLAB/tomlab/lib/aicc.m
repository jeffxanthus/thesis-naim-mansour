% function f = aicc(ssq, N, p)
%
% Compute AICc; The bias corrected Akaike's information criterion (approximated)
%
% ssq   Sum of squares of the residuals
% N     Number of data points
% p     Number of model parameters

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Aug 10, 2002.   Last modified Aug 10, 2002.

function f = aicc(ssq, N, p)

d = N - p - 2;

if d > 0
   f = N*log(ssq/N) + 2*p + 2*(p+1)*(p+2)/d;
else
   % Use standard AIC if too few points
   f = N*log(ssq/N) + 2*p;
end