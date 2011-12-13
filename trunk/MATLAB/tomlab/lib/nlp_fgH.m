% nlp_fgH.m
%
% function [f, g, H] = nlp_fgH(x, Prob, varargin)
%
% nlp_fgH calls the TOMLAB gateway function nlp_f,
% which computes the function value f(x) for the test problem P (Prob.P).
%
% if nargout > 1, nlp_fgH also calls the TOMLAB gateway function nlp_g,
% which computes the gradient g at x, g(x).
%
% if nargout > 2, nlp_fgH also calls the TOMLAB gateway function nlp_H,
% which computes the Hessian H at x, H(x).
%
% nlp_fgH is used when calling OPT TB 2.0.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 6, 1999.    Last modified Sep 9, 1999.

function [f, g, H] = nlp_fgH(x, Prob, varargin)

f=nlp_f(x, Prob, varargin{:});

if nargout > 1
   g=nlp_g(x, Prob, varargin{:});
   if nargout > 2
      H=nlp_H(x, Prob, varargin{:});
   end
end

% MODIFICATION LOG:
%