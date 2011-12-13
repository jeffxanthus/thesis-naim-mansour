% simQG_gdc - gradient and nonlinear constraint vector for simulation
% problem.
%
% function [g, dc] = simQG_gdc(x, Prob)

function [g, dc] = simQG_gdc(x, Prob)

g = [2*(x(1)-x(2));-2*(x(1)-x(2))+3*(x(2)-x(3))^2;-3*(x(2)-x(3))^2+4*(x(3)-x(4))^3;
    -4*(x(3)-x(4))^3+4*(x(4)-x(5))^3;-4*(x(4)-x(5))^3];

dc = [1 0 x(5);2*x(2) 1 0; 3*x(3)^2 -2*x(3) 0;0 1 0;0 0 x(1)]';