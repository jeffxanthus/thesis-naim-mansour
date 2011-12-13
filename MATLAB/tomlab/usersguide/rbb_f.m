% rbb_f - function value for Rosenbrocks Banana, Problem RB BANANA
%
% function f = rbb_f(x, Prob)

function f = rbb_f(x, Prob)

alpha = Prob.user.alpha;

f = alpha*(x(2)-x(1)^2)^2 + (1-x(1))^2;