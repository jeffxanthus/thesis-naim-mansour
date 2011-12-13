% glc5_f.m Compute f(x) for HS 332
%
% function f = glc5_f(x, Prob)

function f = glc5_f(x, Prob)

% HS 332

t = Prob.user.t;
f = pi/3.6*sum((log(t)+x(2)*sin(t)+x(1)*cos(t)).^2+(log(t)+ ...
    x(2)*cos(t)-x(1)*sin(t)).^2);