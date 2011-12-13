% gkls_f.m
%
% function f=gkls_f(x, Prob)
%
% gkls_f evaluates the objective function for test problem P=Prob.P
% at the point x.
%
% The TOMLAB GKLS implementation is based on the GKLS problem generator 
% by M. GAVIANO, D.E. KVASOV, D. LERA, and Ya.D. SERGEYEV.
% http://si.deis.unical.it/~yaro/GKLS.html

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.9.3$
% Written Oct 5, 2005.    Last modified Oct 6, 2005.

function f = gkls_f(x,Prob)

switch(Prob.gklsType)
   case 1,
      % ND
      f = Tgkls('fnd',x);
   case 2,
      % D
      f = Tgkls('fd',x);
   case 3,
      % D2
      f = Tgkls('fd2',x);
end

% MODIFICATION LOG:
%
% 051005 ango Wrote file
% 051006 ango Change to switch statement