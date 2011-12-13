% qpconQG_d2c - The second part of the Hessian to the Lagrangian function for the
%               nonlinear constraints for qpcon quick guide, i.e.:
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x) = H(x) - lam' * d2c(x)
% 
% function d2c = qpconQG_d2c(x, lam, Prob)

function d2c = qpconQG_d2c(x, lam, Prob)

d2c = diag(2*lam(1)./(1+([1:length(x)]-1)/3));