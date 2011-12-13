% function H = ldp_H(x, Prob, varargin)
%
% Computes the Hessian function for the least distance problem
%
%          0.5 * x' * x
%
% The Hessian is I

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written July 17, 2009.  Last modified July 17, 2009.

function H = ldp_H(x, Prob, varargin)

n = length(x);
H = speye(n,n);

% MODIFICATION LOG
%
% 090717  hkh  Written