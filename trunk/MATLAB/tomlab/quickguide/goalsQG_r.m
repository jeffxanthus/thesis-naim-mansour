% function r = goalsQG_r(x, Prob)
%
% Computes residuals to Multi Criterium Unconstrains & Constrained Nonlinear Problem 
% in point x.
% EASY-FIT 'TP355'

function r = goalsQG_r(x, Prob)

r=[11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4); x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)-x(1))];