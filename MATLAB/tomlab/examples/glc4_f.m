% glc4_f: Nonlinear function for Floudas-Pardalos 3.3
%
% function f = glc4_f(x, Prob)
%
% Prob is not used

function f = glc4_f(x, Prob)

f = -25*(x(1)-2)^2-(x(2)-2)^2-(x(3)-1)^2-(x(4)-4)^2-(x(5)-1)^2-(x(6)-4)^2;