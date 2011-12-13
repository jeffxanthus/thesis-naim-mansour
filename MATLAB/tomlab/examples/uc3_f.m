% uc3_f - function value for Generalized Rosenbrocks Banana, 
%         Problem Generalized RB Banana
%
% function f = uc3_f(x, Prob)

function f = uc3_f(x, Prob)

% Safeguarded check on the optional parameter alpha giving steepness in RB
% Default 100, if not given

if isfield(Prob,'userParam')
   if isfield(Prob.userParam,'alpha')
      alpha = Prob.userParam.alpha;
   else
      alpha = 100;
   end
else
   alpha = 100;
end

f = alpha*(x(2)-x(1)^2)^2 + (1-x(1))^2;