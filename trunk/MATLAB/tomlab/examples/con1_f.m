% function f = con1_f(x, Prob)
%
% con1_f computes the objective function f(x)

function f = con1_f(x, Prob)

% Use the standard TOMLAB global variable to be used for communication
% between the function, gradient and Hessian routines

global US_A

% The value could either be directly set as US_A = exp(x(1)); 
% or as structure field:

US_A.e1 = exp(x(1));

% If many items are to be globally sent, structure fields are easier to use.

f = US_A.e1 * (4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1); 