% function c = maximin_c(x,Prob)
%
% Computes constraints of maximin problem: 
%
% max delta s/t ||x-x_i|| >= delta, x_L <= x <= x_U, x in R^n
%
% rewritten as ( delta = x(n+1) )
%
% min -x(n+1) 
% s/t ||x-x_i||_2^2-x(n+1)^2 >= 0, 
%     x_L(1:n) <= x(1:n) <= x_U(1:n), x(n+1) >= 0

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Apr 17, 2005.   Last modified Apr 17, 2005.

function c = maximin_c(x,Prob)

n = length(x)-1;
c = (sum((x(1:n)*ones(1,Prob.mNonLin)-Prob.user.X).^2))'-x(end)^2;
