% jogo09_c.m
%
% function [cx]=jogo09_c(x, Prob)
%
% jogo09_c evaluates the constraints for test problem Prob.P
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 9, 2009.   Last modified June 10, 2009.

function cx = jogo09_c(x, Prob)

P=Prob.P;

switch P
    case 13 % Fully constrained with boundary solutions 5
       cx = [3*(x(1)-2)^2+4*(x(2)-3)^2+2*x(3)^2-7*x(4)-120;
             5*x(1)^2+8*x(2)+(x(3)-6)^2-2*x(4)-40;
             x(1)^2+2*(x(2)-2)^2-2*x(1)*x(2)+14*x(5)-6*x(6);
             0.5*(x(1)-8)^2+2*(x(2)-2)^2+3*x(5)^2-6*x(6)-30;
             -3*x(1)+6*x(2)+12*(x(9)-8)^2-7*x(10)];
    case 14 % Fully constrained with boundary solutions 6
       cx = sum(x.^2);
   otherwise
       cx = [];
end

% MODIFICATION LOG
%
% 090609  bjo  First version created