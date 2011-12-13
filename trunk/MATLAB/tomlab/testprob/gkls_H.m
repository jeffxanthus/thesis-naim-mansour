% gkls_H.m
%
% function H=gkls_H(x, Prob)
%
% gkls_fc evaluates the Hessian of the objective function for test 
% problem P=Prob.P at the point x.
%
% The TOMLAB GKLS implementation is based on the GKLS problem generator 
% by M. GAVIANO, D.E. KVASOV, D. LERA, and Ya.D. SERGEYEV.
% http://si.deis.unical.it/~yaro/GKLS.html

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.6.0$
% Written Oct 5, 2005.    Last modified Oct 18, 2006.

function H = gkls_H(x,Prob)

if Prob.gklsType==3
   H = Tgkls('hd2',x);
else
   error('Hessian not available for ND or D1 GKLS function');
end

% MODIFICATION LOG:
%
% 051005 ango Wrote file
% 061018 ango Error when Hessian not available