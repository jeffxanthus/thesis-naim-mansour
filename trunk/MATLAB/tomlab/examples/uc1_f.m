% uc1_f - function value for Rosenbrocks Banana, Problem RB Banana
%
% function f = uc1_f(x, Prob)

function f = uc1_f(x, Prob)

f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;