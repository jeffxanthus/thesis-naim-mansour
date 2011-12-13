% uc1_g - gradient vector for Rosenbrocks Banana, Problem RB Banana
%
% function g = uc1_g(x, Prob)

function g = uc1_g(x, Prob)

g = [-4*100*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 2*100*(x(2)-x(1)^2)];