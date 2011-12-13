% cpt_c.m
%
% function [cx]=cpt_c(x, Prob)
%
% cpt_c is called by CONOPT to evaluate the nonlinear constraints
% at the point x.
%
% Automatically handles the extended constraint vector when
% both upper and lower bounded constraints are present.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written June 2, 2003.    Last modified May 6, 2004.

function c = cpt_c(x,Prob)

c = feval('nlp_c',x,Prob);
c = c(:);

cexidx = Prob.cexidx;
m2     = Prob.m2;

c = [ c ; c(cexidx(m2+1:end)) ];

% MODIFICATION LOG:
%
% 040506  hkh  Safeguard c, changing to a column vector