% function g = ldp_g(x, Prob, varargin)
%
% Computes the gradient function for the least distance problem
%
%          0.5 * x' * x
%
% The gradient is x


% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written July 17, 2009.  Last modified July 17, 2009.

function g = ldp_g(x, Prob, varargin)

g = x;

% MODIFICATION LOG
%
% 090717  hkh  Written