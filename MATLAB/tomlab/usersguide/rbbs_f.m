% rbbs_f - function value for Constrained Rosenbrocks Banana
%
% function f = rbbs_f(x)

function f = rbbs_f(x)

f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
