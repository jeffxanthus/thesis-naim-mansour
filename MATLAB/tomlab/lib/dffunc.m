% dffunc.m
%
% function [r(j) | c(j)]=dffunc(j, x, Prob)
%
% dffunc calls the TOMLAB gateway functions nlp_r and nlp_c,
% to evaluate the objective and constraint function.
% One function or constraint is evaluated at the time.
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)
%
% dffunc is used as callback from Schittkowski routines

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Jan 18, 2004.   Last modified Jan 18, 2004.

function [rg] = dffunc(j, x, Prob)

global ksNLP_g ksNLP_r ksNLP_x

if isempty(ksNLP_x)
    [ksNLP_r, ksNLP_g] = nlresid(x, Prob);
elseif ~all(ksNLP_x == x)
    [ksNLP_r, ksNLP_g] = nlresid(x, Prob);
end
ksNLP_x = x;

if j > 0
    rg = ksNLP_r(j);
else
    rg = ksNLP_g(-j);
end

% MODIFICATION LOG:
%
% 040118  med  Written, based on nlresid