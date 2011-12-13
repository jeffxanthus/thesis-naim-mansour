% function g = con1_g(x, Prob)
%
% con1_g evaluates the gradient g(x)

function g = con1_g(x, Prob)

global US_A

g = US_A.e1 * [4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+8*x(1)+6*x(2)+1;4*x(1)+4*x(2)+2];