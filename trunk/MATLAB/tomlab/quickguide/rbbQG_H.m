% rbb_H - Hessian matrix for Rosenbrocks Banana, Problem RB BANANA
%
% function H = rbbQH_H(x, Prob)

function H = rbbQG_H(x, Prob)

if isempty(Prob.uP)
   alpha = 100;
else
   alpha = Prob.uP(1);
end

H = [ 12*alpha*x(1)^2-4*alpha*x(2)+2 , -4*alpha*x(1);
                -4*alpha*x(1)        ,    2*alpha  ];