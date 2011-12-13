% cls_r.m
%
% function r = cls_r(x, Prob)
%
% Computes residuals to Constrained Nonlinear Least Squares Problem
% in the point x for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function r = cls_r(x, Prob)

P = Prob.P;
y = Prob.LS.y;
m = size(y,1);

if P==1
    % 'DIST'
    r = x;
elseif P==2
    % 'BAZA'
    r = x-1;
elseif P==3
    % 'GMW'
    r = [x(1); sqrt(2)*x(2)];
elseif P==4
    % 'PAW'
    r = [x(1)-2.4; x(2)-2.1];
elseif P==5
    % 'DAS 1'
    r = [0.1*x(1); x(2)];
elseif P==6
    % 'DAS 2'
    r = [0.75*x(3)+2/3*x(4);
    sqrt(11)/6*x(1)-3/sqrt(11);
    (x(2)-3)/sqrt(2);
    (x(3)+0.5/0.0775)*sqrt(0.0775);
    (x(4)/3-3)/sqrt(2);
    -5/6*x(1)+0.6*x(3)];    
elseif P==7
    % 'DAS 3'
    r = [sqrt(47/136)*(x(1)-4*136/47);
    -x(2)/3+9;
    -x(3)/3+6;
    3/sqrt(17)*(x(1)+17/9*x(2));
    1.5/sqrt(2)*(x(1)+8/9*x(3))];
elseif P==8
    % 'DAS 4'
    r = [sqrt(5)*(x(1)/3-21/5);
    sqrt(7)*(x(2)/4-32/7);
    x(3)-10;
    2*(x(4)-5);
    x(5)-3;
    sqrt(2)*(x(6)-1);
    sqrt(5)*x(7);
    sqrt(7)*(x(8)-11);
    sqrt(2)*(x(9)-10);
    x(10)-7;
    2/3*x(1)+0.75*x(2)];
elseif P==9
    % 'Bob8'
    r = [x(1);
    x(2);
    x(3)];
elseif P==10
    % 'Bob9'
    r = [x(1);
    x(2);
    x(3)];
elseif P==11
    % 'TP001'
    r = [10*(x(2)-x(1)*x(1));
    1-x(1)];
elseif P==12
    % 'TP002'
    r = [10*(x(2)-x(1)*x(1));
    1-x(1)];
elseif P==13
    % 'TP028'
    r = [x(1) + x(2);
    x(2) + x(3)];
elseif P==14
    % 'TP032'
    r=[x(1)+3*x(2)+x(3);
    2*(x(1)-x(2))];
elseif P==15
    % 'TP048'
    r = [x(1)-1;
    x(2)-x(3);
    x(4)-x(5)];
elseif P==16
    % 'TP049'
    r = [x(1)-x(2);
    x(3)-1;
    (x(4)-1)^2;
    (x(5)-1)^3];
elseif P==17
    % 'TP050'
    r = [x(1)-x(2);
    x(2)-x(3);
    (x(3)-x(4))^2;
    x(4)-x(5)];
elseif P==18
    % 'TP051'
    r = [x(1)-x(2);
    x(2)+x(3)-2;
    x(4)-1;
    x(5)-1];
elseif P==19
    % 'TP052'
    r = [4*x(1)-x(2);
    x(2)+x(3)-2;
    x(4)-1;
    x(5)-1];
elseif P==20
    % 'TP053'
    r = [x(1)-x(2);
    x(2)+x(3)-2;
    x(4)-1;
    x(5)-1];
elseif P==21
    % 'TP224'
    r = [sqrt(2)*( x(1)-12 );
    x(2)-20];
elseif P==22
    % 'TP231'
    r = [10*(x(2)-x(1)*x(1));
    1-x(1)];
elseif P==23
    % 'TP269'
    r=zeros(4,1);
    r(1) = x(1)-x(2);
    r(2) = x(2)+x(3)-2;
    r(3) = x(4)-1;
    r(4) = x(5)-1;
elseif P==24
    % 'TP354'
    r=zeros(4,1);
    r(1) = x(1) + 10*x(2);
    r(2) = sqrt(5)*( x(3) - x(4) );
    r(3) = ( x(2)-2*x(3) )^2;
    r(4) = sqrt(10)*( x(1) - x(4) )^2;
elseif (P==25)|(P==26)
    % 'WrHo1 and WrHo2';
    t=Prob.LS.t;
    r = x(1)*exp(-t.*x(5)) + x(2)*exp(-(t - x(9)).^2*x(6)) ...
        + x(3)*exp(-(t - x(10)).^2*x(7)) + x(4)*exp(-(t - x(11)).^2*x(8));
elseif P==27
    % 'RELN'
    r = x-4;
elseif P==28
    % Constrained Walsh
    uP = Prob.uP;
    t = Prob.LS.t;
    r=zeros(m,1);
    if any(x==0)
        r=Inf*ones(m,1);
    else
        for i=1:m,    r(i)=(1-x(1)*t(i)/max(1E-10,x(2)))^(1/(x(1)*uP(1))-1); end
    end
elseif P==29
    %'EASY-TP14'
    r = [x(1)-2;x(2)-1];
elseif P==30
    %'EASY-TP6'
    r = 1-x(1);
elseif P==31
    %'EASY-TP43'
    r = sqrt((x(1)^2+x(2)^2+2*x(3)^2+x(4)^2-5*x(1)-5*x(2)-21*x(3)+7*x(4))+1000);
elseif P==32
    %'EASY-TP57'
    t = Prob.LS.t;
    r = zeros(44,1);
    for i=1:44
        r(i) = x(1)+(0.49-x(1))*exp(-x(2)*(t(i)-8));
    end
elseif P==33
    %'EASY-TP327'
    r = zeros(44,1);
    t = Prob.LS.t;
    for i=1:44
        r(i) = x(1)+(0.49-x(1))*exp(-x(2)*(t(i)-8));
    end
elseif P==34
    %'EASY-TP355'
    r = [11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4); x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)-x(1))];
elseif P==35
    %'EASY-GEO_PROB'
    r = 100-(x(1)^2+x(2)^2+x(3)^2);
elseif P==36
    %'EASY-PSS'
    t = Prob.LS.t;
    if(x(3)==(x(1)+x(2)))
        r = x(4);
    else
        r = zeros(5,1);
        for i=1:5
            r(i) = (x(4)-x(4)/((x(3)-x(1)-x(2))*((x(3)-x(1))*exp(-(x(3)-x(1))*t(i))-x(2)*exp(-x(2)*t(i)))+x(5)));
        end
    end
elseif P==37
    %'EASY-4BAR_LNK'
    t = Prob.LS.t;
    n = 24;
    r1 = zeros(24,1);
    r2 = zeros(24,1);
    alpha1 = zeros(24,1);
    alpha2 = zeros(24,1);
    theta = zeros(24,1);
    pc = zeros(24,1);
    ps = zeros(24,1);
    d = zeros(24,1);
    for i=1:24
        theta(i) = 2*pi*t(i)/n;
        pc(i) = x(1)*cos(theta(i));
        ps(i) = x(1)*sin(theta(i));
        d(i) = x(1)^2+1+2*pc(i);
        beta = pi/4;
        theta_max = 7*pi/9;
        theta_min = pi/6;
        if (x(3)==0 | (x(4)==0 & x(1)==0)) %Avoid divide by zero
            alpha1(i)=1;
            alpha2(i)=1e9;
        else
            alpha1(i) = (x(3)^2+d(i)-x(2)^2)/(2*x(3)*sqrt(d(i))+eps);
            alpha2(i) = ps(i)/(x(4)+pc(i));
        end
        r1(i) = acos(alpha1(i))+atan(alpha2(i));
        r2(i) = 0.5*(theta_max+theta_min)+0.5*(theta_max-theta_min)*sin(theta(i)-beta);
    end
    r = r1-r2;
elseif P==38
    %'EASY-POL_APP
    r = zeros(19,1);
    t = Prob.LS.t;
    for i=1:19
        if(t(i) <= 90)
            r(i) = x(1)+x(2)*t(i)+x(3)*t(i)^2+x(4)*t(i)^3+x(5)*t(i)^4+x(6)*t(i)^5+x(7)*t(i)^6;
        else
            r(i) = x(8)+x(9)*(t(i)-90)+x(10)*(t(i)-90)^2+x(11)*(t(i)-90)^3+x(12)*(t(i)-90)^4+x(13)*(t(i)-90)^5+x(14)*(t(i)-90)^6;
        end
    end
elseif P==39
    %'EASY-LIN_CMP1'
    t = Prob.LS.t;
    r = zeros(11,1);
    if (x(2) ~= x(3))
        a2 = -(x(5) + x(2))/(x(2)-x(3));
    else
        a2 = 0;
    end
    a1 = -1 - a2;
    for i=1:11
        if (t(i) >= x(7))
            r(i) = x(1)*(a1*exp(x(2)*(t(i)-x(7)))+a2*exp(x(3)*(t(i)-x(7))) + 1);
        else
            r = zeros(11,1);
        end
    end
elseif P==40
    %'EASY-DNS'
    t = Prob.LS.t;
    r = zeros(30,1);
    for i=1:30
        r(i) = x(1)*x(2)*(exp(-x(3)*t(i))-exp(-x(2)*t(i)))/(x(2)-x(3));
    end
elseif P==41
    %'EASY-RTD'
    t = Prob.LS.t;
    r = zeros(26,1);
    for i=1:26
        r(i) = 1/(x(2)-x(1))*(exp(-t(i)/x(2))-exp(-t(i)/x(1)));
    end
elseif P==42
    %'EASY-RAT_APP'
    t = Prob.LS.t;
    r = zeros(11,1);
    for i=1:11
        r(i) = x(1)*(t(i)^2 + x(2)*t(i))/(t(i)^2 + x(3)*t(i) + x(4));
    end
elseif P==43
    %'EASY-MAX_PAT1'
    t = Prob.LS.t;
    r = zeros(38,1);
    for i=1:38
        if (x(1) ~= x(2))
            r(i) = 1/(x(2)-x(1))*(exp(-t(i)/x(2))-exp(-t(i)/x(1)));
        else
            r = zeros(38,1);
        end
    end
elseif P==44
    %'EASY-MAX_PAT2'
    t = Prob.LS.t;
    r = zeros(34,1);
    for i=1:34
        if (x(2) ~= x(3) & x(3) ~= 0 & x(2) ~= 0)
            r(i)= x(1)/(x(3)-x(2))*(exp(-t(i)/x(3))-exp(-t(i)/x(2)))+(1-x(1))/x(3)*exp(-t(i)/x(3));
        else
            r = zeros(34,1);
        end
    end
elseif P==45
    %'EASY-TREND'
    a = 971;
    t = Prob.LS.t;
    r = zeros(500,1);
    for i=1:500
        r(i) = x(1)-x(2)*(a-t(i))^x(5)-x(3)*(a-t(i))^x(5)*cos(x(6)*log10(a-t(i)))+x(4)*(a-t(i))^x(5)*sin(x(6)*log10(a-t(i)));
    end
end

if Prob.LS.yUse & m==length(r),  r=r-y; end

% MODIFICATION LOG:
%
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981118  hkh  Use new flag Prob.NLLS.UseYt
% 981210  hkh  Add constrained Walsh problem. Enables convergence with worse x_0
% 001104  hkh  Use field LS, and y, not Yt
% 001105  hkh  Change UseYt to yUse
% 020110  hkh  Add safe guard on Walsh problem (P==28) because of MINOS
% 040517  med  Added 18 problems
% 051128  bjo  Residual for P==42 corrected
% 080603  med  Switched to clsAssign, cleaned