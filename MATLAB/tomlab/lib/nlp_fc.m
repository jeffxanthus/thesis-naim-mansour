% nlp_fc.m
%
% function [f, c]=nlp_fc(x, Prob, varargin)
%
% nlp_fc calls the TOMLAB gateway function nlp_f,
% which computes the function value f(x) for the test problem P (prob.P).
%
% nlp_fc also calls the interface routine nlp_cF, which calls
% the TOMLAB gateway function
% to compute the constraints c(x) for the test problem P (Prob.P).
%
% nlp_fc changes sign on the constraints, giving c(x) <= 0, becuase
% MathWorks Optimization TB is using this formulation.
%
% nlp_fc is used when calling routines in Optimization toolbox (constr.m).
%
% nlp_fc must also explicitely add linear constraints

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 3, 1995.   Last modified Dec 4, 2001.

function [f, c]=nlp_fc(x, Prob, varargin)

f= nlp_f(x, Prob, varargin{:});

[c,Pdim]= nlp_cF( x, Prob, varargin{:});

if isempty(c)
   % Bug in constr, have to return something for c
   c=0;
else
   c=-c;
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 011204  hkh  Set c=0 to avoid crash when no constraints
%
% In Pdim the different constraint dimensions are stored
%
% (1) Total number of constraints (not bounds) (=mcon)
% (2) Number of linear and nonlinear equalities (=me)
% (3) Number of linear equalities (=meq)
% (4) Number of linear inequalities (=mleq)