% nsdp_d2c.m:
%
% function d2c = nsdp_d2c(x_k, lam, Prob, varargin)
%
% TOMLAB gateway routine for computation of function values d2c(x) in
% nonlinear semidefinite programming.
%
% nsdp_d2c calls the user routine in Prob.USER.d2c with the standard variables
% and matrix variables separated as vector x and cell array of matrices Y.
%
% Bjorn Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written October 8, 2008. Last modified October 8, 2008.

function d2c = nsdp_d2c(x, lam, Prob, varargin)

% Update x-values in matrices Y{i}, i=1:k
[x, Y] = nsdp_x2xY(x,Prob.PENOPT.NSDP);

d2c = feval(Prob.USER.d2c, x, Y, lam, Prob, varargin);

% MODIFICATION LOG
%
% 081008 bjo  Wrote file