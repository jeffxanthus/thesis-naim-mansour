% rbb_g - gradient vector for Rosenbrocks Banana, Problem RB BANANA
%
% function g = rbb_g(x, Prob)

function g = rbb_g(x, Prob)

alpha = Prob.user.alpha;

g = [-4*alpha*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 2*alpha*(x(2)-x(1)^2)];