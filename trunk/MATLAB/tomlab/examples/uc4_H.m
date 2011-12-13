% uc4_H - Hessian value for simple problem x_1^2 + (x_2-2)^4
%
% function H = uc4_H(x, Prob)

function H = uc4_H(x, Prob)

H = [ 2 0; 0,  12*(x(2)-2)^2  ];