% function dc = maximin_dc(x,Prob)
%
% Computes constraint Jacobians of maximin problem: 
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

function dc = maximin_dc(x,Prob)

m         = Prob.mNonLin;
n         = Prob.N-1;
dc        = zeros(m,n+1);
dc(:,end) = -2*x(end);
dc(:,1:n) = 2*(ones(m,1)*x(1:n)'-Prob.user.X');
