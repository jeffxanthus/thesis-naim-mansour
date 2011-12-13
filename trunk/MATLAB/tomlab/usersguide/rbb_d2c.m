% rbb_d2c - The second part of the Hessian to the Lagrangian function for the
%     nonlinear constraints for Rosenbrocks Banana, Problem RB BANANA, i.e.
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x) = H(x) - lam' * d2c(x)
% 
% function d2c=crbb_d2c(x, lam, Prob)

function d2c=rbb_d2c(x, lam, Prob)

% The only nonzero element in the second derivative matrix for the single
% constraint is the (1,1) element, which is a constant -2.

d2c = lam(1)*[-2 0; 0 0];

