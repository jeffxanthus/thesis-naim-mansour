% con_H.m
%
% function H=con_H(x, Prob)
%
% con_H computes the Hessian to the objective function f at the point x for
% test problem Prob.P.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function H=con_H(x, Prob)

global Baux

P=Prob.P;

if (P==1)|(P==2)|(P==3)
   % Exponential problem 1,2,3
   e1=exp(x(1));
   a=e1*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+16*x(1)+10*x(2)+9);
   b=e1*(4*x(2)+4*x(1)+6);
   H=[a b;b 4*e1];
elseif P==4
   % Powell 1969
   H = zeros(5,5);
   for i=1:5
       for j=1:5
           prod=1;
           prodi=1;
           prodj=1;
           prodij=1;
           for k=1:5
               prod=prod*x(k);
               if i~=k
                  prodi=prodi*x(k);
               end
               if j~=k
                  prodj=prodj*x(k);
               end
               if (i~=k)&(j~=k)
                  prodij=prodij*x(k);
               end
           end
           if i==j
              H(i,j)=0;
           else
              H(i,j)=prodij*exp(prod);
           end
           H(i,j)=H(i,j)+prodi*prodj*exp(prod);
       end
   end
elseif P==5
   % Fletcher 12.19
   H=zeros(2,2);
elseif P==6
   % Circle-triangle
   H=zeros(4,4);
   H(1,2)=1;
   H(2,1)=1;
elseif P==7
   % Chvatal
   N=Prob.uP(1);
   if N <= 0
      N=5;
   end
   H=zeros(N,N);
elseif P==8
   % Schittkowski 14
   H=2*eye(2);
elseif P==9
   % Schittkowski 24
   r=9*sqrt(3);
   a=2*x(2)^3/(3*r);
   b=2*(x(1)-3)*x(2)^2/r;
   H=[a b; b 2*x(2)*((x(1)-3)^2-9)/r];
elseif P==10
   % Schittkowski 66
   H=zeros(3,3);
elseif P==11
   % Fletcher page 279. Penalty algorithm has slow convergence.
   H=zeros(2,2);
elseif P==12
   % Nonlinear maximation for ABB Robotics, project spring 1996
   H=[2*Baux,zeros(6); zeros(6,12)];
elseif P==13
   % Hock-Shittkowski 375
   H = -2*eye(10);
elseif P==14
   % DAS 2
   if isfield(Prob.PartSep,'index')
      i=Prob.PartSep.index;
   else
      i=0;
      Prob.PartSep.pSepFunc=0;
   end
   if i==1 & Prob.PartSep.pSepFunc==6 % Compute only 1st term
      H = zeros(4,4);
      H(1,1) = 11/36;
   elseif i==2 & Prob.PartSep.pSepFunc==6 % Compute only 2nd term
      H = zeros(4,4);
      H(2,2) = 1/2;
   elseif i==3 & Prob.PartSep.pSepFunc==6 % Compute only 3rd term
      H = zeros(4,4);
      H(3,3) = 0.0775;
   elseif i==4 & Prob.PartSep.pSepFunc==6 % Compute only 4th term
      H = zeros(4,4);
      H(4,4) = 1/18;
   elseif i==5 & Prob.PartSep.pSepFunc==6 % Compute only 5th term
      H = [   25/36    0     -0.6*5/6    0     
                0      0         0       0
            -0.6*5/6   0       0.36      0
                0      0         0       0];
   elseif i==6 & Prob.PartSep.pSepFunc==6 % Compute only 6th term
      H = [     0      0         0       0     
                0      0         0       0
                0      0       0.5625   0.5
                0      0        0.5     4/9];
   else
      H = [    1        0    -0.6*5/6    0
               0       0.5       0       0
           -0.6*5/6     0        1      0.5 
               0        0       0.5     0.5];
   end
elseif P==15
   % Entropy problem
   H = sparse(diag(1 ./ max(1E-300,x)));
elseif P==16
   H = [ 2-1/x(1)^2  0 ; 0 2-1/x(2)^2 ];
elseif P==17
   H = zeros(4,4);
end

% MODIFICATION LOG
%
% 030101  hkh  Error in P=9, unsymmetric Hessian
% 030123  hkh  Added entropy problem
% 030128  hkh  Safe guard x for entropy problem, avoid division with 0
% 041102  ango Added problem 16
% 050104  ango Added problem 17 - BMI as NLP
% 080603  med  Switched to conAssign, cleaned