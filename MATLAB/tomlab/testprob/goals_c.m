% function cx=goals_c(x, Prob)
%
% The goals problems are 
% Multi Criterium Unconstrained & Constrained Nonlinear Problems
%
% goals_c evaluates the nonlinear constraints for test problem P (Prob.P) 
% at the point x.
 
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function cx=goals_c(x, Prob)

P=Prob.P;
cx=[];

if P==1
   %'EASY-FIT TP269'
elseif P==2
   %'EASY-TP13'
   cx = (1-x(1))^3-x(2);
elseif P==3
   %'EASY-FIT TP46'
   cx = [x(1)^2*x(4)+sin(x(4)-x(5)); x(2)+((x(3)^4)*x(4)^2)];
elseif P==4
   %'EASY-FIT TP48'
elseif P==5
   %'EASY-FIT TP354'
elseif P==6
   %'EASY-FIT TP355'
   c1 = 11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4);
   c2 = x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)- x(1));
   c3 = 11-4*x(1)*x(4)-4*x(2)*x(4)+x(3)*x(4);
   c4 = 2*x(1)+20*x(2)-0.5*x(3)+2*x(4)+2*x(2)*x(4)*(x(3)-4*x(1));
   cx = (c1)^2+(c2)^2-(c3)^2-(c4)^2;
elseif P==7
   %'EASY-FIT TP372'
   cx = zeros(12,1);
   cx(1)  = x(1)+x(2)*exp(-5*x(3))+x(4);
   cx(2)  = x(1)+x(2)*exp(-3*x(3))+x(5);
   cx(3)  = x(1)+x(2)*exp(-x(3))+x(6);
   cx(4)  = x(1)+x(2)*exp(x(3))+x(7);
   cx(5)  = x(1)+x(2)*exp(3*x(3))+x(8);
   cx(6)  = x(1)+x(2)*exp(5*x(3))+x(9);
   cx(7)  = -x(1)-x(2)*exp(-5*x(3))+x(4);
   cx(8)  = -x(1)-x(2)*exp(-3*x(3))+x(5);
   cx(9)  = -x(1)-x(2)*exp(-x(3))+x(6);
   cx(10) = -x(1)-x(2)*exp(x(3))+x(7);
   cx(11) = -x(1)-x(2)*exp(3*x(3))+x(8);
   cx(12) = -x(1)-x(2)*exp(5*x(3))+x(9);
elseif P==8
   %'EASY-FIT TP373'
   cx = zeros(6,1);
   cx(1) = x(1)+x(2)*exp(-5*x(3))+x(4);
   cx(2) = x(1)+x(2)*exp(-3*x(3))+x(5);
   cx(3) = x(1)+x(2)*exp(-x(3))+x(6);
   cx(4) = x(1)+x(2)*exp(x(3))+x(7);
   cx(5) = x(1)+x(2)*exp(3*x(3))+x(8);
   cx(6) = x(1)+x(2)*exp(5*x(3))+x(9);
end

% MODIFICATION LOG:
%
% 040517  med  Created.
% 050503  hkh  Vectorized 7,8, clean up
% 080603  med  Switched to clsAssign, cleaned