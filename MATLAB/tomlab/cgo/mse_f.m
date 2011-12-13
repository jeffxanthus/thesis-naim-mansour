% function mse = mse_f(x, Prob)
%
% Kriging variance, (less complicated variant), Sasena Thesis pg. 58
%
% mse(y^(x)) = sigma^2(x) = 
%     sigma^2_z * (1 - r_x' R^(-1) r_x  +  (1-1'R^(-1) r_x)^2/(1'R^(-1)1) )
%                                                        
%
% sigma^2_z is a constant, Prob.EGO.sigma2
% 
% To be maximized, change sign

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written May 2, 1998.   Last modified May 3, 2005.

function mse = mse_f(x, Prob)

invR   = Prob.EGO.invR;
k      = Prob.EGO.k;
p      = Prob.EGO.p;
n      = size(invR,1);

if length(p) > 1
   % p in R^d, same length as x, all p in ]0,2]
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^(p*ones(1,n))) ))';
else
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^p) ))';
end

%r = zeros(n,1);
%for i = 1:n
%   r(i) = exp(-(  theta'*(abs( x-X(:,i) )).^p ));
%end

Rr = invR * r;

% MSE. Minus sign to make it a minimization problem
mse = -Prob.EGO.sigma2*(1-r'*Rr + (1-sum(Rr))^2/sum(sum(invR)));

% MODIFICATIONS
% 050502  hkh  Written, with code from ego_c
% 050503  hkh  Expanded with additional term, similar to EGO paper
