% minlpQG_d2c - The second part of the Hessian to the Lagrangian function for the
%     nonlinear constraints for minlp quick guide, i.e.:
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x) = H(x) - lam' * d2c(x)
% 
% function d2c=minlpQG_d2c(x, lam, Prob)

function d2c=minlpQG_d2c(x, lam, Prob)

d2c      = spalloc(5,5,2);
d2c(1,1) = lam(1)*2;
d2c(2,2) = lam(2)*3/(4*sqrt(x(2)));