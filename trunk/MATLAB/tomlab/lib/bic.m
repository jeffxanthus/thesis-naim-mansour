% function f = bic(ssq, N, p)
%
% Compute BIC; The Bayesian information criterion (approximated)
%
% ssq   Sum of squares of the residuals
% N     Number of data points
% p     Number of model parameters

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Aug 10, 2002.   Last modified Aug 02, 2005.

function f = bic(ssq, N, p)

f = N*log(ssq/N) + p + p*log(N);

% MODIFICATION LOG
%
% 050801  med  Function name changed to bic