% dfgrad.m
%
% function [Jdg]=dfgrad(j, x, Prob)
%
% dfgrad calls the TOMLAB gateway functions nlp_J and nlp_dc,
% to evaluate the Jacobian of the objective and constraint function.
% One objective and constraint is evaluated at the time.
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)
%
% dfgrad is used as callback from Schittkowski routines

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Jan 18, 2004.   Last modified Apr 7, 2004.

function [Jdg] = dfgrad(j, x, Prob)

global ksNLP_dg ksNLP_J ksNLP_xJ

if isempty(ksNLP_xJ)
    [ksNLP_J, ksNLP_dg]=nlJac(x, Prob);
elseif ~all(ksNLP_xJ == x)
    [ksNLP_J, ksNLP_dg]=nlJac(x, Prob);
end
ksNLP_xJ = x;

if j > 0
    Jdg = ksNLP_J(j,:);
else
    Jdg = ksNLP_dg(-j,:);
end

% MODIFICATION LOG:
%
% 040118  med  Written, based on dffunc
% 040407  ango Return entire rows