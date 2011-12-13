% uc1_H - Hessian matrix for Rosenbrocks Banana, Problem RB Banana
%
% function H = uc1_H(x, Prob)

function H = uc1_H(x, Prob)

H = [ 1200*x(1)^2-400*x(2)+2 , -400*x(1);
                 -400*x(1)   ,  200  ];