% nlp_gdc.m
%
% function [g, dc]=nlp_gdc(x, Prob, varargin)
%
% nlp_gdc calls the TOMLAB gateway function nlp_g,
% which computes the gradient of f at x for the test problem P (prob.P).
%
% nlp_gdc also calls nlp_dcF, which in turn calls the TOMLAB gateway function
% nlp_dc, that computes the gradient for all constraints at x,
% dc, for the test problem P.
%
% nlp_gdc changes sign on the constraints, giving c(x) <= 0, because
% MathWorks Optimization TB is using this formulation.
%
% nlp_gdc also explicitely add the derivative of the linear constraints
%
% nlp_gdc is used when calling routines in Optimization toolbox (constr.m).
%
% Number of rows should be equal to the number of variables

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 3, 1995.   Last modified Dec 4, 2001.

function [g, dc]=nlp_gdc(x, Prob, varargin)

nargin;
g= nlp_g(x, Prob, varargin{:});
dc=-full(nlp_dcF(x, Prob, varargin{:}));

if isempty(dc)
    dc = ones(length(x),1);
else
    dc = dc';
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990628  hkh  Add init to call for nlp_dcF
% 990705  hkh  Make dc a full matrix
% 011204  hkh  Transpose dc and create dummy dc for unconstrained problems