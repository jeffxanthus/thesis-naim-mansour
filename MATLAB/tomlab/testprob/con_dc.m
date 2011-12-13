% con_dc.m
%
% function dc=con_dc(x, Prob)
%
% con_dc computes the gradient to the constraints c at the point x for
% the test problem Prob.P
%
% One row for each constraint, one column for each variable

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function dc=con_dc(x, Prob)

global B1 B2 B3 B4 B5 B6 mass_matrix

P=Prob.P;

if P==1
   % Exponential problem 1
   dc=[-x(2)+1 -x(1)+1; x(2) x(1)];
elseif P==2
   % Exponential problem 2
   dc=[-x(2)+1 -x(1)+1; x(2) x(1)];
elseif P==3
   % Exponential problem 3
   dc=[-x(2)+1 -x(1)+1; x(2) x(1)];
elseif P==4
   %  Powell 1969
   dc=zeros(3,length(x));
   for i = 1:length(x)
       dc(1,i)=2*x(i);
   end
   dc(2,2:5)=[x(3), x(2),-5*x(5),-5*x(4)];
   dc(3,1:2)=[3*x(1)^2,3*x(2)^2];
elseif P==5
   %  Fletcher 12.19
   dc=[-2*x(1),1];
elseif P==6
   % Circle-triangle
   h1=x(1)^2+x(2)^2;
   h2=x(1)*x(3)+x(2)*x(4);
   if h1==0 
      dc=[1E10,1E10,1E10,1E10];
   else
      dc=[2*h2*x(3)/h1-2*x(1)*(h2/h1)^2, 2*h2*x(4)/h1-2*x(2)*(h2/h1)^2,...
          2*h2*x(1)/h1-2*x(3), 2*h2*x(2)/h1-2*x(4)];
   end
elseif P==7
   N=Prob.uP(1);
   if N <= 0
      N=5;
   end
   dc=-eye(N,N);
   for i=1:N
       for j=i+1:N
           dc(i,j)=-2*10^(j-i);
       end
   end
elseif P==8
   %  Schittkowski 14
   dc=[-0.5*x(1), -2*x(2)];
elseif P==9
   % Schittkowski 24
   dc=[];
elseif P==10
   % Schittkowski 66
   dc=[-exp(x(1)) 1 0; 0 -exp(x(2)) 1];
elseif P==11
%  Fletcher page 279. Penalty algorithm has slow convergence.
   dc=[2*x(1),2*x(2)];
elseif P==12
   % Nonlinear maximation for ABB Robotics, VT 1996
   dc=zeros(length(x)/2,length(x));
   dc(1,1:6) = 2*x(1:6)'*B1;
   dc(2,1:6) = 2*x(1:6)'*B2;
   dc(3,1:6) = 2*x(1:6)'*B3;
   dc(4,1:6) = 2*x(1:6)'*B4;
   dc(5,1:6) = 2*x(1:6)'*B5;
   dc(6,1:6) = 2*x(1:6)'*B6;
   dc(:,7:12) = mass_matrix;
elseif P==13
   % Hock-Shittkowski 375
   dc = 2*x'./(1+([1:length(x)]-1)/3);
elseif P==14
   % DAS 2
   dc = [];
elseif P==15
   % Entropy
   dc = [];
elseif P==16
   dc = [ x(2) x(1) ];
elseif P==17
   dc = [2-2*x(4)  6     0     -2*x(1) ; ...
         0        -4  8-2*x(4) -2*x(3) ];
end

% MODIFICATION LOG:
% 041102  ango Added problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 080603  med  Switched to conAssign, cleaned