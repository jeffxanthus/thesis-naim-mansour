% gkls_g.m
%
% function g=gkls_g(x, Prob)
%
% gkls_g evaluates the gradient of the objective function for test 
% problem P=Prob.P at the point x.
%
% The TOMLAB GKLS implementation is based on the GKLS problem generator 
% by M. GAVIANO, D.E. KVASOV, D. LERA, and Ya.D. SERGEYEV.
% http://si.deis.unical.it/~yaro/GKLS.html

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.6.0$
% Written Oct 5, 2005.    Last modified Oct 18, 2006.

function g = gkls_g(x,Prob)

switch(Prob.gklsType)
   case 1,
      error('Gradient not available for ND GKLS function'); 
   case 2,
      g = Tgkls('gd',x);
   case 3,
      g = Tgkls('gd2',x);
end

% MODIFICATION LOG:
%
% 051005 ango Wrote file
% 051006 ango Change to switch statement
% 061018 ango Error when gradient not available