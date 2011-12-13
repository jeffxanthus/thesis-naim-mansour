% function d2c=con_d2c(x, lam, Prob)
%
% con_d2c computes the 2nd part of the Hessian to the Lagrangian function,
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function d2c=con_d2c(x, lam, Prob)

global B1 B2 B3 B4 B5 B6

P=Prob.P;
uP=Prob.uP;

if P==1
   % Exponential problem 1
   d2c=[         0  [-1 1]*lam; 
        [-1 1]*lam           0]; 
elseif P==2
   % Exponential problem 2
   d2c=[         0  [-1 1]*lam; 
        [-1 1]*lam           0]; 
elseif P==3
   % Exponential problem 3
   d2c=[           0  [-1 1]*lam; 
        [-1 1]*lam             0]; 
elseif P==4
   %  Powell 1969
   d2c=diag(lam(1)*[2 2 2 2 2]);
   d2c(2,3)=lam(2); d2c(4,5)=-5*lam(2);
   d2c(3,2)=lam(2); d2c(5,4)=-5*lam(2);
   d2c(1,1)=d2c(1,1)+6*x(1)*lam(3);
   d2c(2,2)=d2c(2,2)+6*x(2)*lam(3);
elseif P==5
   %  Fletcher 12.19
   d2c=[-2*lam(1) 0; 0 0]; 
elseif P==6
   % Circle-triangle
   if lam(1) ~=0
      d2c=Inf*ones(4,4);
      if x(1)==0 & x(2)==0
         return
      end
      h1=x(1)^2+x(2)^2;
      h2=x(1)*x(3)+x(2)*x(4);
      d2c(1,1)=2*x(3)^2/h1-4*x(1)*x(3)*h2/h1^2-...
           2*(h2/h1)^2-4*x(1)*x(3)*h2/h1^2+8*x(1)^2*h2^2/h1^3;
      d2c(1,2)=2*x(3)*x(4)/h1-4*x(2)*x(3)*h2/h1^2-...
           2*(h2/h1)^2-4*x(1)*x(4)*h2/h1^2-8*x(1)*x(2)*h2^2/h1^3;
      d2c(1,3)=2*x(1)*x(3)/h1+2*h2/h1-4*x(1)^2*h2/h1^2;
      d2c(1,4)=2*x(2)*x(3)/h1-4*x(1)*x(2)*h2/h1^2;
      d2c(2,2)=2*x(4)^2/h1-4*x(2)*x(4)*h2/h1^2-...
           2*(h2/h1)^2-4*x(2)*x(4)*h2/h1^2+8*x(2)^2*h2^2/h1^3;
      d2c(2,3)=2*x(1)*x(4)/h1-4*x(1)*x(2)*h2/h1^2;
      d2c(2,4)=2*x(2)*x(4)/h1+2*h2/h1-4*x(2)^2*h2/h1^2;
      d2c(3,3)=2*x(1)^2/h1-2;
      d2c(4,4)=2*x(2)^2/h1-2;
      d2c(3,4)=2*x(1)*x(2)/h1;
      d2c(2,1)=d2c(1,2);
      d2c(3,1)=d2c(1,3);
      d2c(4,1)=d2c(1,4);
      d2c(3,2)=d2c(2,3);
      d2c(4,2)=d2c(2,4);
      d2c(4,3)=d2c(3,4);
      d2c=lam(1)*d2c;
   else
      d2c=zeros(4,4);
   end
elseif P==7
   %  Chvatal
   if uP(1) ==0
      N=5;
   else
      N=uP(1);
   end
   d2c=zeros(N,N);
elseif P==8
   %  Schittkowski 14
   d2c=zeros(2,2);
   if lam(1) > 0
      d2c(1,1)=-0.5*lam(1);
      d2c(2,2)=-2*lam(1);
   end
elseif P==9
   % Schittkowski 24
   d2c=zeros(2,2);
elseif P==10
   % Schittkowski 66
   d2c=zeros(3,3);
   d2c(1,1)=-exp(x(1))*lam(1);
   d2c(2,2)=-exp(x(2))*lam(2);
elseif P==11
%  Fletcher page 279. Penalty algorithm has slow convergence.
   d2c=2*lam(1)*eye(2);
elseif P==12
   % Nonlinear maximation for ABB Robotics, VT 1996
   d2c=zeros(12,12);
   d2c(1:6,1:6)=2*(B1*lam(1)+B2*lam(2)+B3*lam(3)+B4*lam(4)+B5*lam(5)+B6*lam(6));
elseif P==13
   % Hock-Shittkowski 375
   d2c = diag(2*lam(1)./(1+(1:length(x)-1)/3));
elseif P==14
   % DAS 2
   d2c = [];
elseif P==15
   % Entropy
   d2c = [];
elseif P==16
   d2c = lam(1)*[ 0 1 ; 1 0 ];
elseif P==17
   d2c = spalloc(4,4,2);
   d2c(1,4) = -2*lam(1);
   d2c(3,4) = -2*lam(2);
   d2c = d2c+d2c';
end

% MODIFICATION LOG:
%
% 981116  hkh  Error in problem 3, extra zeros.
% 981120  hkh  Error in problem 9, lam in exponent.
% 990831  hkh  Fixed circle-triangle example
% 030113  ango Error in problem 4, diagonal entries were overwritten
% 041102  ango Added problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 050610  ango Fixed empty output for P==15 (was missing)
% 080603  med  Switched to conAssign, cleaned