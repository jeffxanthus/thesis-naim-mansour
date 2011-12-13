% simQG_fc - function value and constraints for simulation problem.
%
% function [f, c] = simQG_fc(x, Prob)

function [f, c] = simQG_fc(x, Prob)

f = (x(1)-x(2))^2+(x(2)-x(3))^3+(x(3)-x(4))^4+(x(4)-x(5))^4;

c = [x(1)+x(2)^2+x(3)^3-3;x(2)-x(3)^2+x(4)-1;x(1)*x(5)-1];