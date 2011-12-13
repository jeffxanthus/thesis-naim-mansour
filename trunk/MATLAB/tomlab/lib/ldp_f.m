% function f = ldp_f(x, Prob, varargin)
%
% Computes the objective function for the least distance problem
%
%          0.5 * x' * x

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written July 17, 2009.  Last modified July 17, 2009.

function f = ldp_f(x, Prob, varargin)

f = 0.5*(x'*x);

% MODIFICATION LOG
%
% 090717  hkh  Written