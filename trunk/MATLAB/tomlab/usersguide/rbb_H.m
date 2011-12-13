% rbb_H - Hessian matrix for Rosenbrocks Banana, Problem RB BANANA
%
% function H = crbb_H(x, Prob)

function H = rbb_H(x, Prob)

alpha = Prob.user.alpha;

H = [ 12*alpha*x(1)^2-4*alpha*x(2)+2 , -4*alpha*x(1);
                -4*alpha*x(1)        ,    2*alpha  ];