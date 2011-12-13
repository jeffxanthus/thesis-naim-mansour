% function H = uc_H(x,Prob)
%
% Hessian matrix for test functions for Unconstrained optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Aug 4, 2008.

function H = uc_H(x,Prob)

x=x(:);
P=Prob.P;
uP=Prob.uP;

if P == 1
    if isfield(Prob.PartSep,'index')
        i=Prob.PartSep.index;
    else
        i=0;
        Prob.PartSep.pSepFunc=0;
    end
    if i==1 & Prob.PartSep.pSepFunc==2
        H = [ 1200*x(1)^2-400*x(2),-400*x(1);
            -400*x(1),            200];
    elseif i==2 & Prob.PartSep.pSepFunc==2
        H = [ 2 0; 0 0];
    else
        H = [ 1200*x(1)^2-400*x(2)+2,-400*x(1);
            -400*x(1),              200];
    end
elseif P == 2
    e=exp(x(1));
    e_H11 = e*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 16*x(1) + 10*x(2) + 9);
    e_H12 =e*(4*x(1) + 4*x(2) + 6);
    H = [e_H11,e_H12;
        e_H12,4*e];
elseif P == 3
    v = x(1) * x(2) * x(3);
    if real(v)~=0
        e = 200 / v;
        H = [  2*e/x(1)^2      50+ e/(x(1)*x(2))  100+ e/(x(1)*x(3))
            50+e/(x(1)*x(2)) 2*e/x(2)^2         200+ e/(x(2)*x(3))
            100+e/(x(1)*x(3)) 200+ e/(x(2)*x(3)) 2*e/x(3)^2 ];
    else
        H = [  400      50  100
            50     400  200
            100     200  400 ];
    end
elseif P == 4
    r = norm(real(x));
    if r==0, r=1E-10; end
    if real(x(2)) >= 0
        phi = acos(x(1)/r);
    else
        phi = 2*pi - acos(x(1)/r);
    end
    c=x(1)/r;
    s=x(2)/r;
    sr = sin(phi + r);
    cr = cos(phi + r);
    a = c - abs(s)/r;
    if real(s)==0
        b = 0;
        h12 =  cr * (-c*b)/r;
    else
        b = s*(1 + c/(abs(s)*r));
        h12 = -a*sr*b + cr * (s*abs(s)/r-c*(b+2*c*s/(r*abs(s))))/r;
    end
    H2 = uP(2)*[ s^2/r -s*c/r; -s*c/r c^2/r];
    H = uP(1)*[-a^2*sr + cr * s*(b-3*c/r)/r h12;...
        h12  -b^2*sr + cr * c*(a-abs(s)/r)/r] + H2;
elseif P == 5
    H = 12*[x(1)-x(2)-0.5,-x(1)+x(2)+0.5; -x(1)+x(2)+0.5,x(1)];
elseif P == 6
    H = [2 0 ; 0 8];
elseif P == 7	% AG-1, Anders Goeran 1.
    % 2 min at f(0.354667,-0.354667)=f(-0.354667,0.354667)=1.42052
    e1=exp(x(1)^2+x(2)^2);
    e2=exp(2*x(1)*x(2)-1);
    h12=4*x(1)*x(2)*e1+(2+4*x(1)*x(2))*e2+2;
    H = [(2+4*x(1)^2)*e1+4*x(2)^2*e2  h12;h12 (2+4*x(2)^2)*e1+4*x(1)^2*e2];
elseif P == 8	% AG-2, Anders Goeran 2.
    % min at f(0,0)=-1. Deep hole, flat function away from origo.
    e1=exp(-x(1)^2-x(2)^2);
    h12=-4*x(1)*x(2)*e1;
    H = [(2-4*x(1)^2)*e1  h12;h12 (2-4*x(2)^2)*e1];
elseif P == 9% RE-1
    %Marsden - Tromba fig. 2.1.17, min at f(0,0)=0, max at(0,+-1)=3.
    eVal=exp(1-x(1)^2-x(2)^2);
    e1= uP(1)*eVal;   % uP(1) is sign
    h11=2*(1-5*x(1)^2-3*x(2)^2+2*x(1)^4+6*(x(1)*x(2))^2);
    h21=-4*x(1)*x(2)*(4-x(1)^2-3*x(2)^2);
    h22=2*(3-x(1)^2-15*x(2)^2+6*x(2)^4+2*(x(1)*x(2))^2);
    H = e1*[h11 h21;h21 h22];
elseif P == 10% RE-2
    %netlib/trig, n=2, several max, min and saddle points.
    H = [cos(x(1)) 0; 0 3*cos(x(2))+sin(x(2))];
elseif P == 11% Structured exponential, only total hessian provided
    eX=exp(x(1));
    if uP(1) ==0
        H11 = 4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 16*x(1) + 10*x(2) + 9;
        H12 = 4*x(1) + 4*x(2) + 6;
        H = eX*[H11,H12; H12, 4];
    elseif uP(1)==1
        H = [eX*(4*x(1)^2 + 16*x(1) + 8) 0; 0 0 ];
    elseif uP(1)==2
        H = eX*[2*x(2)^2 4*x(2); 4*x(2) 4];
    elseif uP(1)==3
        H = eX*[4*x(1)*x(2) + 8*x(2), 4*x(1) + 4; 4*x(1) + 4, 0];
    elseif uP(1)==4
        H = eX*[2*x(2) + 1, 2; 2 0];
    end
elseif P == 12	% Fletcher Q.2.4
    H = [ 12*x(1)^2+20*x(1)^3-4*x(2) -4*x(1) ; -4*x(1) 2];
elseif P == 13	% Fletcher Q.2.5
    H=[ 4+12*x(1)+12*x(1)^2 -2; -2 2 ];
elseif P == 14	% Fletcher Q.2.2
    % min i f(0.69,-1.37)=-0.58
    H = [12*x(1)^2  1;1  2];
elseif P == 15	% Fletcher Q.2.7 Quadratic
    H = [2  0;0  8];
elseif P == 16	% Fletcher Q.2.6
    H = [ 2*x(2)^2-8*x(2)+8 ...
        4*x(1)*x(2)-8*x(1)+4*x(2)-8;0 2*x(1)^2+4*x(1)+2];
    H(2,1)=H(1,2);
elseif P == 17	% Fletcher Q.3.3
    e = exp(x(1)^2-x(2)^2);
    H = [1+5*x(1)^2+2*x(1)^2*x(2)^2+x(2)^2+2*x(1)^4, ...
        -2*x(1)*x(2)*(x(1)^2+x(2)^2);
        0 , 1-x(1)^2+2*x(1)^2*x(2)^2-5*x(2)^2+2*x(2)^4 ];
    H(2,1)=H(1,2);
    H = e*H;
elseif P == 18  % Nash-Sofer page 307 #3
    H=[60*x(1)^2-12 2;2 72*x(2)^2+10];
elseif P == 19  % 3rd degree 
    H=[10-12*x(1) -2*x(2);-2*x(2) 2-2*x(1)];
end

% MODIFICATION LOG:
%
% 981021  hkh  Added handling of Partial Separable functions to P==1
% 981102  hkh  Check if pSepFunc == 2
% 080603  med  Switched to conAssign, cleaned
