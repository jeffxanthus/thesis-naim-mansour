% function [r,J] = lsq_rJ(x, Prob)
%
% Example of low level function to compute residual and Jacobian to be used
% when running lsqnonlin
%
% Three problems from testprob\ls_prob.m (#2,#4,#5)
%
% The problems illustrate how to send specific problem information to the 
% ser functions
% In the third problem (Population problem) also (t,y) are defined before
% optimization and sent to this routine.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 28, 1999.  Last modified July 28, 1999.

function [r,J] = lsq_rJ(x, Prob)

if nargin < 2
   P=1;
   uP=[];
else
   P  = Prob.P;
   uP = Prob.uP;
end

x=x(:);

if P==1
   % Walsh
   % C = uP(1) 
   if isempty(uP)
      uP=96.05;
   end
   t=[ 2000; 5000; 10000; 20000; 30000; 50000];
   y=[0.9427; 0.8616; 0.7384; 0.5362; 0.3739; 0.3096];
   m=size(t,1);
   r=zeros(m,1);
   if any(x==0)
      r=Inf*ones(m,1);
   else
      for i=1:m,    r(i)=(1-x(1)*t(i)/x(2))^(1/(x(1)*uP(1))-1); end
   end
   r=r-y;
elseif P==2
   % Gisela
   % K = uP(1)
   if isempty(uP)
      uP=5;
   end
   t=[0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
      96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];
   y=[30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
      19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];
   m=size(t,1);
   r=zeros(m,1);
   for i=1:m
      r(i)=uP(1)*x(1)*(exp(-x(2)*t(i))-exp(-x(1)*t(i)))/(x(3)*(x(1)-x(2)));
   end
   r=r-y;
elseif P==3
   % Population problem
   
   % Pick up t and y from Prob structure
   t=Prob.LS.t;
   y=Prob.LS.y;
   m=size(t,1);
   r=zeros(m,1);
   for i=1:m,    r(i)=x(1)*x(2)^t(i); end
   r=r-y;
end

if nargout > 1
  % Compute Jacobian
  if P==1
     % Walsh
     if x(1)==0
        b=1E100;
     else
        b=1/(x(1)*uP(1))-1;
     end
     J=zeros(m,2);
     for i=1:m
         a=1-x(1)*t(i)/x(2);
         if a <= 0
            i
            t(i)
            a
         end
         if a <= 0, loga=1E100; else loga=log(a);end
         if x(2)==0
            J(i,1)=1E100;
            J(i,2)=1E100;
         elseif x(1)==0
            J(i,1)=1E100;
            J(i,2)=a^(b-1)*(x(1)*t(i)/x(2)^2)*b;
         else
            J(i,1)=a^b*(-1/(uP(1)*x(1)^2)*loga-b*t(i)/(x(2)*a));
            J(i,2)=a^(b-1)*(x(1)*t(i)/x(2)^2)*b;
         end
     end
  elseif P==2
     % Gisela: uP(1) = K
     a=uP(1)*x(1)/(x(3)*(x(1)-x(2)));
     b=x(1)-x(2);
     J=zeros(m,3);
     for i=1:m
         e1=exp(-x(1)*t(i)); e2=exp(-x(2)*t(i));
         J(i,1)=a*(t(i)*e1+(e2-e1)*(1-1/b));
         J(i,2)=a*(-t(i)*e2+(e2-e1)/b);
         J(i,3)=-a*(e2-e1)/x(3);
     end
  elseif P==3
     % Population problem
     J=zeros(m,2);
     for i=1:m, J(i,:)=[x(2)^t(i),t(i)*x(1)*x(2)^(t(i)-1)]; end
  end
end

% MODIFICATION LOG:
%

