% con_f.m
%
% function f=con_f(x, Prob)
%
% con_f computes the objective function f in the point x
% for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function f=con_f(x, Prob)

global US_A mass_vector Baux gravity

P=Prob.P;

if (P==1) | (P==2) | (P==3)
   % Exponential problem 1,2,3
   f=exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
elseif P==4
   % Powell 1969
   f=exp(x(1)*x(2)*x(3)*x(4)*x(5));
elseif P==5
   % Fletcher 12.19
   f=x(1)+x(2);
elseif P==6
   % Circle-Triangle
   f=x(1)*x(2);
elseif P==7
   % Chvatal
   N=Prob.uP(1);
   if N <= 0
      N=5;
   end
   f=0;
   for j=1:N
       f=f-x(j)*10^(j-1);
   end
elseif P==8
   % Schittkowski 14
   f=(x(1)-2)^2+(x(2)-1)^2;
elseif P==9
   % Schittkowski 24
   f=((x(1)-3)^2-9)*x(2)^3/(27*sqrt(3));
elseif P==10
   % Schittkowski 66
   f=0.2*x(3)-0.8*x(1);
elseif P==11
   % Fletcher page 279. Penalty algorithm has slow convergence.
   f=-x(1)-x(2);
elseif P==12
   % Nonlinear maximation for ABB Robotics, project spring 1996
   %    x(1:6) is the angle speed of each axis
   %    x(7:12) is the angle acceleration of each axis 
   %    The speed is limited to -beta <= x(1:6) <= beta
   %
   % nonlinear constraints (quadratic):
   %    -alfa_i <= tau_i <= alfa_i
   %
   % In con_prob.m the global variables are defined, e.g. defining
   % mass_vector, Baux, gravity
   f=x(1:6)'*Baux*x(1:6) + mass_vector*x(7:12) + gravity;
elseif P==13
   % Hock-Shittkowski 375
   f = -x'*x;
elseif P==14
   % DAS 2
   if isfield(Prob.PartSep,'index')
      i=Prob.PartSep.index;
   else
      i=0;
      Prob.PartSep.pSepFunc=0;
   end
   if i==1 & Prob.PartSep.pSepFunc==6 % Compute only 1st term
      f = 0.5*( sqrt(11)/6*x(1)-3/sqrt(11) )^2; 
   elseif i==2 & Prob.PartSep.pSepFunc==6 % Compute only 2nd term
      f = 0.5*( ( x(2)-3)/sqrt(2) )^2;
   elseif i==3 & Prob.PartSep.pSepFunc==6 % Compute only 3rd term
      f = 0.5*0.0775*(x(3)+0.5/0.0775)^2;
   elseif i==4 & Prob.PartSep.pSepFunc==6 % Compute only 4th term
      f = 0.5*( (x(4)/3-3)/sqrt(2) )^2;
   elseif i==5 & Prob.PartSep.pSepFunc==6 % Compute only 5th term
      f = 0.5*( -5/6*x(1)+0.6*x(3) )^2;
   elseif i==6 & Prob.PartSep.pSepFunc==6 % Compute only 6th term
      f = 0.5*( 0.75*x(3)+2/3*x(4) )^2;
   else
      f = 0.5*( sqrt(11)/6*x(1)-3/sqrt(11) )^2 + ...
          0.5*( (x(2)-3)/sqrt(2) )^2 + ...
          0.5*( (x(3)+0.5/0.0775)*sqrt(0.0775) )^2 + ...
          0.5*( (x(4)/3-3)/sqrt(2) )^2 + ...
          0.5*( -5/6*x(1)+0.6*x(3) )^2 + ...
          0.5*( 0.75*x(3)+2/3*x(4) )^2;
   end
elseif P==15
   % Entropy problem
   US_A = log(max(1E-300,x));
   f    = sum( x.*US_A );
elseif P==16
   f = x(1)^2 + x(2)^2 + log(x(1)*x(2));
elseif P==17
   f = -x(4);
end

% MODIFICATION LOG
%
% 030123  hkh  Added entropy problem
% 030128  hkh  Safe guard x in entropy problem avoiding log(0)
% 041102  ango Add SOCS 6.3 example, Problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 080603  med  Switched to conAssign, cleaned