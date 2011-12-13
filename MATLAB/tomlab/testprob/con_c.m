% con_c.m
%
% function cx=con_c(x, Prob)
%
% con_c evaluates the constraints for test problem P=Prob.P
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function cx=con_c(x, Prob)

global mass_matrix B1 B2 B3 B4 B5 B6 gravity_c

P=Prob.P;

if P==1
   % Exponential problem 1
   cx=[-1.5-x(1)*x(2)+x(1)+x(2); 10+x(1)*x(2)];
elseif P==2
   % Exponential problem 2. Added lower bound zero for both parameters
   cx=[-1.5-x(1)*x(2)+x(1)+x(2); 10+x(1)*x(2)];
elseif P==3
   % Exponential problem 3. Added equality, but no lower bounds zero
   cx=[-1.5-x(1)*x(2)+x(1)+x(2); 10+x(1)*x(2)];
elseif P==4
   % Powell 1969
   cx=[x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2; x(2)*x(3)-5*x(4)*x(5);...
       x(1)^3+x(2)^3];
elseif P==5
   % Fletcher 12.19
   cx=x(2)-x(1)^2;
elseif P==6
   % Circle-Triangle
   if x(1)==0 & x(2)==0
      cx=1;
   else
      cx=(x(1)*x(3)+x(2)*x(4))^2/(x(1)^2+x(2)^2)-x(3)^2-x(4)^2;
   end
elseif P==7
   % Chvatal. Lower bound set as x_L
   N=Prob.uP(1);
   if N <= 0
      N=5;
   end
   cx=zeros(N,1);
   for i=1:N
       s=0;
       for j=i+1:N
           s=s+10^(j-i)*x(j);
       end
       cx(i)=10^(2*N-2*i)-x(i)-2*s;
   end;
elseif P==8
   % Schittkowski 14. Bracken, McCormick, Himmelblau. Start (2,2), f=1
   % Min in (.5(sqrt(7)-1), .25(sqrt(7)+1))
   % f(x^*)=9 - 2.875*sqrt(7)
   cx=-0.25*x(1)^2-x(2)^2;
elseif P==9
   % Schittkowski 24. Betts and Box. Start in (1,.5), f= -.01336459
   % Min in (3, sqrt(3)). f(x^*)=-1
   cx=zeros(0,1);
elseif P==10
   % Schittkowski 66. Eckhardt. Start in (0,1.05,2.9). f=.58
   % Min in (.1841264879, 1.202167873,3.327322322), f(x^*)=.5181632741
   cx=[x(2)-exp(x(1)); x(3)-exp(x(2))];
elseif P==11
   % Fletcher page 279. Penalty algorithm has slow convergence.
   cx=x(1)^2+x(2)^2;
elseif P==12
   % Nonlinear maximation for ABB Robotics, VT 1996
   cx(1) = mass_matrix(1,:)*x(7:12)+x(1:6)'*B1*x(1:6) + gravity_c(1);
   cx(2) = mass_matrix(2,:)*x(7:12)+x(1:6)'*B2*x(1:6) + gravity_c(2);
   cx(3) = mass_matrix(3,:)*x(7:12)+x(1:6)'*B3*x(1:6) + gravity_c(3);
   cx(4) = mass_matrix(4,:)*x(7:12)+x(1:6)'*B4*x(1:6) + gravity_c(4);
   cx(5) = mass_matrix(5,:)*x(7:12)+x(1:6)'*B5*x(1:6) + gravity_c(5);
   cx(6) = mass_matrix(6,:)*x(7:12)+x(1:6)'*B6*x(1:6) + gravity_c(6);
   cx=cx(:);
elseif P==13
   % Hock-Shittkowski 375
   cx = sum(x.^2./(1+([1:length(x)]'-1)/3));
elseif P==14
   % DAS 2
   cx = [];
elseif P==15
   % Entropy
   cx = [];
elseif P==16
   cx = x(1)*x(2);
elseif P==17
   cx = [ 2*x(1)+6*x(2)-2*x(1)*x(4) ; ...
         -4*x(2)+8*x(3)-2*x(3)*x(4) ];
end

% MODIFICATION LOG
%
% 981018  hkh  Found error in Circle/Triangle when x=(0,0,*,*)
% 041102  ango Added problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 080603  med  Switched to conAssign, cleaned