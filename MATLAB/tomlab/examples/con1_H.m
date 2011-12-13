% function H = con1_H(x, Prob)
%
% con1_H evaluates the gradient H(x)

function H = con1_H(x, Prob)

global US_A

a = 4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 16*x(1) + 10*x(2) + 9;
b = 4*x(2) + 4*x(1) + 6;

H = US_A.e1 * [a b;b 4];