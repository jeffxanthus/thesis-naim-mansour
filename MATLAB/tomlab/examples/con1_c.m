% function c = con1_c(x, Prob)
%
% con1_c evaluates the constraints c(x)

function c = con1_c(x, Prob)

% Two nonlinear (quadratic) constraints

c = [ - x(1)*x(2) + x(1) + x(2); x(1)*x(2)];