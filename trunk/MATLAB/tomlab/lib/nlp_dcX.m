% nlp_dcX.m
%
% Calls the interface routine nlp_dcF
% which computes the constraint gradients dc at x for test problem P (Prob.P).
%
% nlp_dcX also explicitely adds linear constraints
%
% function [dc]=nlp_dcX(x, Prob, varargin)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Aug 25, 1998.  Last modified Oct 19, 2000.

function [dc]=nlp_dcX(x, Prob, varargin)

global n_dc

n_dc=n_dc+1;

dc= nlp_dcF( x, Prob, varargin{:});