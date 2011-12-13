% function [f,J,domerr]=mcp_vifunjac(x,jacflag,Prob) 
%
% TOMLAB gateway routine for the computation of PATH information.
%
% INPUT:
% x       Current iterate x.
% jacflag Indicates whether or not the Jacobian is needed.
% Prob    TOMLAB problem structure.
%
% OUTPUT:
% F       Function evaluations at x.
% J       The Jacobian of F.
% domerr  The number of divided by zero (set to 0).
%
% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 18, 2004.  Last modified Apr 18, 2004.

function [f,J,domerr]=mcp_vifunjac(x,jacflag,Prob) 

fJ = Prob.PATH.fJ;
m = Prob.PATH.m;
n = Prob.N;
A = Prob.PATH.A;
b = Prob.PATH.b;

if jacflag == 1
  [f,J,domerr] = feval(fJ,x(1:n),jacflag,Prob);
  if (m > 0)
    f = [f - A'*x(n+1:n+m); 
         A*x(1:n) - b];
    J = [J, -A'; A, sparse(m,m)];
  end
else
  [f,J,domerr] = feval(fJ,x(1:n),jacflag,Prob);
  f = [f - A'*x(n+1:n+m); 
       A*x(1:n) - b];
end