% ls_r.m
%
% function r = ls_r(x, Prob)
%
% Computes residuals to Nonlinear Least Squares problem in the point x 
% for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function r = ls_r(x, Prob)

global LS_A

x=x(:);
P  = Prob.P;
uP = Prob.uP;
y  = Prob.LS.y;
t  = Prob.LS.t;

m=size(y,1);
if P==1
   % Powell
   r=[x(1); uP(1)*x(1)^2+x(1)];
elseif P==2
   % Walsh
   r=zeros(m,size(x,2));
   if any(x==0)
      r=Inf*ones(m,1);
   else
      for i=1:m,    r(i)=(1-x(1)*t(i)/x(2))^(1/(x(1)*uP(1))-1); end
   end
elseif P==3
   % Gisela
   r=zeros(m,size(x,2));
   for i=1:m
      r(i)=uP(1)*x(1)*(exp(-x(2)*t(i))-exp(-x(1)*t(i)))/(x(3)*(x(1)-x(2)));
   end
elseif P==4
   m=round(length(y(:))/2);
   r=zeros(2*m,size(x,2));
   t=zeros(2*m,size(x,2)); % Use special t for this problem
   y=reshape(y(:),m,2);
   for i=1:m
       z=atan((x(3)-y(i,2))/max( (x(2)-y(i,1)),1E-8 ));
       % atan gives solution in -pi/2 to pi/2, correct if outside 0 to pi/2
       if (y(i,2)-x(3))*z < 0
          z=z+pi;
       end
       t(i)=cos(z);
       r(i)=x(1)*t(i)+x(2);
       k=m+i;
       t(k)=sin(z);
       r(k)=x(1)*t(k)+x(3);
   end
   % May send t with global, but otherwise must recompute in ls_J.m
   LS_A=t;
elseif P==5
   % Population problem
   r=zeros(m,size(x,2));
   for i=1:m,    r(i)=x(1)*x(2)^t(i); end
elseif P==6
   % Plasmid n=2
   p_0 = uP(2);
   r=zeros(m,size(x,2));
   for i=1:m
      denom = ( p_0*x(1) + x(2) )*exp( (x(1)+x(2))*t(i) ) - x(2)*(1-p_0);
      numer = ( p_0*x(1) + x(2) )*exp( (x(1)+x(2))*t(i) ) + x(1)*(1-p_0);
      r(i) = denom/numer;
   end
elseif P==7
   % Plasmid n=3
   r=zeros(m,size(x,2));
   for i=1:m
      e1 = exp((x(1)+x(2))*t(i));
      r(i) = ((x(3)*x(1)+x(2))*e1-x(2)*(1-x(3)))/...
             ((x(3)*x(1)+x(2))*e1+x(1)*(1-x(3)));
   end   
elseif P==8
   % Plasmid n=3 (subst.)
   D = uP(3);
   r=zeros(m,size(x,2));
   for i=1:m
      denom = (D-x(2)+x(3))*exp((x(1)-x(2)+x(3))*t(i))-x(3)*...
              (1-(D-x(2))/(x(1)-x(2)));
      numer = (D-x(2)+x(3))*exp((x(1)-x(2)+x(3))*t(i))+x(1)-D;
      r(i) = denom/numer;
   end   
elseif P==9
   % Plasmid n=3 (probability)
   D = uP(3);
   r=zeros(m,size(x,2));
   for i=1:m
      e1 = exp((x(1)+x(2)*(x(3)-1))*t(i));
      d1 = (D+x(2)*(x(3)-1));
      r(i) = 1-(x(3)*x(2)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)/(d1*e1+x(1)-D);
   end   
elseif P==10
   % Parameterized test function (Huschens)
   phi = uP(1); 
   r=zeros(m,size(x,2));
   r(1) = x(1);
   r(2) = (x(1)-2*phi)*x(2);
   r(3) = x(2);
elseif P==11
   % Signomial problem
   n = uP(1); 
   l = Prob.user.l;
   r=-ones(m,size(x,2));
   C = Prob.user.C;
   A = Prob.user.A;
   for k =1:l
       z = ones(m,size(x,2));
       for j=1:n
           z=z.*x(j).^A(:,j,k);
       end
       r = r + C(:,k).*z;
   end
elseif P==12
   % Signomial problem, pseudorand
   n = uP(1); 
   l = Prob.user.l;
   r=-ones(m,size(x,2));
   C = Prob.user.C;
   A = Prob.user.A;
   for k =1:l
       z = C(:,k);
       for j=1:n
           z=z.*x(j).^A(:,j,k);
       end
       r = r + z;
   end
elseif P==13
   % Exponential problem
   r=-ones(m,size(x,2));
   l = Prob.user.l;
   C = Prob.user.C;
   A = Prob.user.A;
   for k =1:l
       r = r + C(:,k).*exp(A(:,:,k)*x);
   end
elseif P==14
   % Exponential problem
   r=-zeros(m,size(x,2));
   l = Prob.user.l;
   C = Prob.user.C;
   A = Prob.user.A;
   for k =1:l
       r = r + C(:,k).*exp(A(:,:,k)*x);
   end
elseif P==15
   % Trigonometric problem
   r = (-Prob.user.C+Prob.user.A*sin(x)+Prob.user.B*cos(x)).^2;
end
if Prob.LS.yUse & m==length(r),  r=r-y; end

% MODIFICATION LOG:
%
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981108  hkh  Use general global variable LS_A instead of circle_t
% 981118  hkh  Use new flag Prob.NLLS.UseYt
% 990314  hkh  Remove unnecessary yMod computation
% 030323  hkh  Adding four problems
% 031201  hkh  Use zeros(m,size(x,2)) to make AD tool MAD work, size(x,2)=1
% 050302  hkh  Added problems 11-15
% 080603  med  Switched to clsAssign, cleaned