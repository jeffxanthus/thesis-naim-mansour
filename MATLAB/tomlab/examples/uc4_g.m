% uc4_g - gradient value for simple problem x_1^2 + (x_2-2)^4
%
% function g = uc4_g(x, Prob)

function g = uc4_g(x, Prob)

g = [2*x(1); 4*(x(2)-2)^3];