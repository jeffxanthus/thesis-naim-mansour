% nlp_JT.m
%
% function JT=nlp_JT(x, Prob, varargin)
%
% Computes the transposed Jacobian for a nonlinear least squares problem
% Used in calls to MathWorks Optimization toolbox, routines leastsq and conls.
%
% Calls the TOMLAB gateway routine nlp_J to compute the Jacobian matrix

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 15, 1995. Last modified Oct 19 2000.

function JT=nlp_JT(x, Prob, varargin)

JT = full(nlp_J(x, Prob, varargin{:})');
if isempty(JT)
   JT = zeros(length(x),0);
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990706  hkh  Make JT full matrix