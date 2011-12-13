% con_g.m
%
% function g=con_g(x, Prob)
%
% con_df computes the gradient to objective function f in the point x for the
% test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function g=con_g(x, Prob)

global US_A mass_vector Baux

P=Prob.P;

if (P==1)|(P==2)|(P==3)
   % Exponential problem 1,2,3
   e1=exp(x(1));
   g=[e1*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+8*x(1)+6*x(2)+1);e1*(4*x(1)+4*x(2)+2)];
elseif P==4
   % Powell 1969
   e5=exp(x(1)*x(2)*x(3)*x(4)*x(5));
   g=[x(2)*x(3)*x(4)*x(5)*e5; x(1)*x(3)*x(4)*x(5)*e5; x(1)*x(2)*x(4)*x(5)*e5;...
      x(1)*x(2)*x(3)*x(5)*e5; x(1)*x(2)*x(3)*x(4)*e5];
elseif P==5
   % Fletcher 12.19
   g=[1;1];
elseif P==6
   % Circle-triangle
   g=[x(2);x(1);0;0];
elseif P==7
   % Chvatal
   N=Prob.uP(1);
   if N <= 0
      N=5;
   end
   g=zeros(N,1);
   for j=1:N
       g(j)=-10^(j-1);
   end
elseif P==8
   % Schittkowski 14
   g=[2*x(1)-4;2*x(2)-2];
elseif P==9
   % Schittkowski 24
   r=9*sqrt(3);
   g=[2*(x(1)-3)*x(2)^3/(3*r);x(2)^2*((x(1)-3)^2-9)/r];
elseif P==10
   % Schittkowski 66
   g=[-0.8 0 0.2]';
elseif P==11
%  Fletcher page 279. Penalty algorithm has slow convergence.
   g=[-1;-1];
elseif P==12
   % Nonlinear maximation for ABB Robotics, project spring 1996
   g = [2*Baux*x(1:6); mass_vector'];
elseif P==13
   % Hock-Shittkowski 375
   g = -2*x;
elseif P==14
   % DAS 2
   if isfield(Prob.PartSep,'index')
      i=Prob.PartSep.index;
   else
      i=0;
      Prob.PartSep.pSepFunc=0;
   end
   if i==1 & Prob.PartSep.pSepFunc==6 % Compute only 1st term
      g = [11/36*x(1)-0.5;0;0;0]; 
   elseif i==2 & Prob.PartSep.pSepFunc==6 % Compute only 2nd term
      g = [0;(x(2)-3)/2;0;0]; 
   elseif i==3 & Prob.PartSep.pSepFunc==6 % Compute only 3rd term
      g = [0;0;0.0775*(x(3)+0.5/0.0775);0]; 
   elseif i==4 & Prob.PartSep.pSepFunc==6 % Compute only 4th term
      g = [0;0;0;(x(4)/3-3)/6]; 
   elseif i==5 & Prob.PartSep.pSepFunc==6 % Compute only 5th term
      g = [-5/6*( -5/6*x(1)+0.6*x(3) );0;0.6*( -5/6*x(1)+0.6*x(3) );0]; 
   elseif i==6 & Prob.PartSep.pSepFunc==6 % Compute only 6th term
      g = [0;0;0.75*(0.75*x(3)+2/3*x(4));2/3*(0.75*x(3)+2/3*x(4))];
   else
      g = [11/36*x(1)-0.5-5/6*(-5/6*x(1)+0.6*x(3));...
          (x(2)-3)/2;...
          0.0775*(x(3)+0.5/0.0775)+0.6*(-5/6*x(1)+0.6*x(3))+0.75*(0.75*x(3)+2/3*x(4));...
          (x(4)/3-3)/6+2/3*(0.75*x(3)+2/3*x(4))];
   end
elseif P==15
   % Entropy problem
   g = 1 + US_A;
elseif P==16
   g = [ 2*x(1) + 1/x(1) ; ...
         2*x(2) + 1/x(2) ];
elseif P==17
   g = [0 0 0 -1]';
end

% MODIFICATION LOG
%
% 030123  hkh  Added entropy problem
% 041102  ango Added problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 080603  med  Switched to conAssign, cleaned