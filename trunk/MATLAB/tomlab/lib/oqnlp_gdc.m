% function [g,dc] = oqnlp_gdc(x, Prob, varargin)
%
% TOMLAB interface routine to avoid numerical differentiation of
% integer variables in MINLP problems
%
% Computation of the gradient vector and the constraint
% gradient matrix at the same time

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 12, 2004.   Last modified Nov 12, 2004.

function [g,dc] = oqnlp_gdc(x, Prob, varargin)

n          = Prob.N;
IntVars    = Prob.MIP.IntVars;
g          = NaN*ones(n,1);
g(IntVars) = 0;
m          = Prob.mNonLin;
if m > 0
    dc            = NaN*Prob.ConsPattern;
    dc(:,IntVars) = 0;
else
    dc            = [];
end

% MODIFICATION LOG
%
% 041112  hkh  Written