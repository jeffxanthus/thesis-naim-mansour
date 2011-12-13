% glc4_c: Nonlinear constraints for Floudas-Pardalos 3.3
%
% function c = glc4_c(x, Prob)

function c = glc4_c(x, Prob)

% Two nonlinear constraints (quadratic)

c = [(x(3)-3)^2+x(4); (x(5)-3)^2+x(6)];