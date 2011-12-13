% rbb_f - function value for Rosenbrocks Banana, Problem RB BANANA
%
% function f = rbbQG_f(x, Prob)

function f = rbbQG_f(x, Prob)

if isempty(Prob.uP)
   alpha = 100;
else
   alpha = Prob.uP(1);
end

f = alpha*(x(2)-x(1)^2)^2 + (1-x(1))^2;