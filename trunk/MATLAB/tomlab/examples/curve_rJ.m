% function [r,J] = curve_rJ(x, t, Prob)
%
% Example of low level function to compute residual and Jacobian to be used
% when running lsqcurvefit
%
% Three problems from ls_prob (#2,#4,#5)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 28, 1999.  Last modified Apr 22, 2005.

function [r,J] = curve_rJ(x, t, Prob)


if nargin < 3
   P=1;
   uP=[];
else
   P  = Prob.P;
   uP = Prob.uP;
end

x=x(:);

m=size(t,1);

if P==1
   % Walsh
   % C = uP(1) 
   if isempty(uP)
      uP=96.05;
   end
   r=zeros(m,1);
   if any(x==0)
      r=Inf*ones(m,1);
   else
      for i=1:m,    r(i)=(1-x(1)*t(i)/x(2))^(1/(x(1)*uP(1))-1); end
   end
elseif P==2
   % Gisela
   % K = uP(1)
   if isempty(uP)
      uP=5;
   end
   r=zeros(m,1);
   for i=1:m
      r(i)=uP(1)*x(1)*(exp(-x(2)*t(i))-exp(-x(1)*t(i)))/(x(3)*(x(1)-x(2)));
   end
elseif P==3
   % Population problem
   r=zeros(m,1);
   for i=1:m,    r(i)=x(1)*x(2)^t(i); end
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
         %if a <= 0
         %   i
         %   t(i)
         %   a
         %end
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
% 050422 hkh Avoid printout of a, i, t(i), when bad x for Walsh

