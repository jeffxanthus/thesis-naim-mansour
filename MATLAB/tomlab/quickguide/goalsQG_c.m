% function cx=goalsQG_c(x, Prob)
%
% The goals problems are Multi Criterium Unconstrained & Constrained Nonlinear Problems
% 'EASY-FIT TP355'
%
% goalsQG_c evaluates the nonlinear constraints at the point x.

function cx=goalsQG_c(x, Prob)

c1 = 11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4);
c2 = x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)- x(1));
c3 = 11-4*x(1)*x(4)-4*x(2)*x(4)+x(3)*x(4);
c4 = 2*x(1)+20*x(2)-0.5*x(3)+2*x(4)+2*x(2)*x(4)*(x(3)-4*x(1));
cx = [(c1)^2+(c2)^2-(c3)^2-(c4)^2];