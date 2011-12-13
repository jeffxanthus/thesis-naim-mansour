% nsdp_g.m:
%
% function g = nsdp_g(x_k, Prob, varargin)
%
% TOMLAB gateway routine for computation of function gradient values g(x) 
% in nonlinear semidefinite programming.
%
% nsdp_g calls the user routine in Prob.USER.g with the standard variables
% and matrix variables separated as vector x and cell array of matrices Y.
%
% Bjorn Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written October 8, 2008. Last modified October 8, 2008.

function g = nsdp_g(x, Prob, varargin)

% Update x-values in matrices Y{i}, i=1:k
[x, Y] = nsdp_x2xY(x,Prob.PENOPT.NSDP);

g = feval(Prob.USER.g, x, Y, Prob, varargin);

% MODIFICATION LOG
%
% 081008 bjo  Wrote file