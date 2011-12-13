% rbb_g - gradient vector for Rosenbrocks Banana, Problem RB BANANA
%
% function g = rbbQG_g(x, Prob)

function g = rbbQG_g(x, Prob)

if isempty(Prob.uP)
   alpha = 100;
else
   alpha = Prob.uP(1);
end

g = [-4*alpha*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 2*alpha*(x(2)-x(1)^2)];