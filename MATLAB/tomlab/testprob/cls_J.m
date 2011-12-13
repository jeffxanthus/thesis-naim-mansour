% cls_J.m
%
% function J = cls_J(x, Prob)
%
% Computes the Jacobian to the Constrained Nonlinear Least Squares Problem
% in the point x for the test problem P (Prob.P).
%
% J(i,j) is dr_i/d_x_j

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function J = cls_J(x, Prob)

P=Prob.P;
m=size(Prob.LS.y,1);

if P==1
    % 'DIST'
    J = eye(2);
elseif P==2
    % 'BAZA'
    J = eye(4);
elseif P==3
    % 'GMW'
    J = [ 1    0
        0 sqrt(2) ];
elseif P==4
    % 'PAW'
    J = eye(2);
elseif P==5
    % 'DAS 1'
    J = [0.1 0; 0 1];
elseif P==6
    % 'DAS 2'
    J = eye(m,4);
    J(1,1)=sqrt(11)/6;
    J(2,2)=1/sqrt(2);
    J(3,3)=sqrt(0.0775);
    J(4,4)=J(2,2)/3;
    J(5,1)=-5/6;
    J(5,2)=0;
    J(5,3)=0.6;
    J(5,4)=0;
    J(6,1)=0;
    J(6,2)=0;
    J(6,3)=0.75;
    J(6,4)=2/3;
elseif P==7
    % 'DAS 3'
    J = zeros(m,3);
    J(1,1) = sqrt(47/136);
    J(2,2) = -1/3;
    J(3,3) = J(2,2);
    J(4,1) = 3/sqrt(17);
    J(4,2) = 1/J(4,1);
    J(4,3) = 0;
    J(5,1) = 1.5/sqrt(2);
    J(5,2) = 0;
    J(5,3) = sqrt(8)/3;
elseif P==8
    % 'DAS 4'
    J = eye(m,10);
    J(1,1) = sqrt(5)/3;
    J(2,2) = sqrt(7)/4;
    J(4,4) = 2;
    J(6,6) = sqrt(2);
    J(7,7) = sqrt(5);
    J(8,8) = sqrt(7);
    J(9,9) = J(6,6);
    J(11,1)= 2/3;
    J(11,2)= 0.75;
elseif P==9
    % 'Bob8'
    J = eye(3);
elseif P==10
    % 'Bob9'
    J = eye(3);
elseif P==11
    % 'TP001'
    J = [-20*x(1)  10
        -1      0 ];
elseif P==12
    % 'TP002'
    J = [-20*x(1)  10
        -1      0 ];
elseif P==13
    % 'TP028'
    J = [1 1 0
        0 1 1 ];
elseif P==14
    % 'TP032'
    J = [1  3  1
        2 -2  0];
elseif P==15
    % 'TP048'
    J = [ 1  0  0  0  0
        0  1 -1  0  0
        0  0  0  1 -1 ];
elseif P==16
    % 'TP049'
    J = [ 1  -1   0     0        0
        0   0   1     0        0
        0   0   0  2*(x(4)-1)  0
        0   0   0     0    3*(x(5)-1)^2 ];
elseif P==17
    % 'TP050'
    J = [ 1   -1        0               0          0
        0    1       -1               0          0
        0    0   2*(x(3)-x(4))  -2*(x(3)-x(4))   0
        0    0        0               1         -1 ];
elseif P==18
    % 'TP051'
    J = [ 1 -1  0  0  0
        0  1  1  0  0
        0  0  0  1  0
        0  0  0  0  1 ];
elseif P==19
    % 'TP052'
    J = [ 4 -1  0  0  0
        0  1  1  0  0
        0  0  0  1  0
        0  0  0  0  1 ];
elseif P==20
    % 'TP053'
    J = [ 1 -1  0  0  0
        0  1  1  0  0
        0  0  0  1  0
        0  0  0  0  1 ];
elseif P==21
    % 'TP224'
    J = [ sqrt(2)  0
        0      1 ];
elseif P==22
    % 'TP231'
    J = [-20*x(1)  10
        -1      0 ];
elseif P==23
    % 'TP269'
    J = [ 1 -1  0  0  0
        0  1  1  0  0
        0  0  0  1  0
        0  0  0  0  1 ];
elseif P==24
    % 'TP354'
    J = [1  10  0  0
        0  0  sqrt(5)  -sqrt(5)
        0  2*(x(2)-2*x(3))  -4*(x(2)-2*x(3))  0
        2*sqrt(10)*(x(1)-x(4))  0  0  -2*sqrt(10)*(x(1)-x(4))];
elseif (P==25)|(P==26)
    % 'WrHo1 and WrHo2';
    J=zeros(m,11);
    t=Prob.LS.t;
    for i = 1 : 65
        t9   = t(i) - x(9);
        t10  = t(i) - x(10);
        t11  = t(i) - x(11);
        s9   = t9 * t9;
        s10  = t10 * t10;
        s11  = t11 * t11;
        e1   = exp(-t(i)*x(5));
        e2   = exp(-s9*x(6));
        e3   = exp(-s10*x(7));
        e4   = exp(-s11*x(8));
        r2   = x(2)*e2;
        r3   = x(3)*e3;
        r4   = x(4)*e4;
        J(i,1) = e1;
        J(i,2) = e2;
        J(i,3) = e3;
        J(i,4) = e4;
        J(i,5) = -t(i)*x(1)*e1;
        J(i,6) = -s9*r2;
        J(i,7) = -s10*r3;
        J(i,8) = -s11*r4;
        J(i,9) =  2*t9*x(6)*r2;
        J(i,10) =  2*t10*x(7)*r3;
        J(i,11) =  2*t11*x(8)*r4;
    end
elseif P==27
    % 'RELN'
    J = eye(m);
elseif P==28
    % Walsh
    % uP(1) = C
    uP = Prob.uP;
    t = Prob.LS.t;
    if x(1)==0
        b=1E100;
    else
        b=1/(x(1)*uP(1))-1;
    end
    J=zeros(m,2);
    for i=1:m
        % Add safe guard for division with 0.
        % Solvers have lower bound x_L(2) = 1E-10, but MINOS calls with 0.
        a=1-x(1)*t(i)/max(1E-10,x(2));
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
elseif P==29
    %'EASY-TP14'
    J = eye(2);
elseif P==30
    %'EASY-TP6'
    J = [-1 0];
elseif P==31
    %'EASY-TP43'
    f = sqrt((x(1)^2+x(2)^2+2*x(3)^2+x(4)^2-5*x(1)-5*x(2)-21*x(3)+7*x(4))+1000);
    J = [(2*x(1)-5)/(2*f) (2*x(2)-5)/(2*f) (4*x(3)-21)/(2*f) (2*x(4)+7)/(2*f)];
elseif P==32
    %'EASY-TP57'
    t=Prob.LS.t;
    J = zeros(44,2);
    for i=1:44
        J(i,:) = [1-exp(-x(2)*(t(i)-8)) (0.49-x(1))*(-t(i)+8)*exp(-x(2)*(t(i)-8))];
    end
elseif P==33
    %'EASY-TP327'
    J = zeros(44,2);
    t = Prob.LS.t;
    for i=1:44
        J(i,:) = [1-exp(-x(2)*(t(i)-8))  -(t(i)-8)*(0.49-x(1))*exp(-x(2)*(t(i)-8))];
    end
elseif P==34
    %'EASY-TP355'
    J = [-x(4)  -x(4)  x(4) -x(1)-x(2)+x(3);
        1-x(2)*x(4)  10+x(4)*x(3)-x(1)*x(4) -1+x(2)*x(4)  1+x(2)*x(3)-x(2)*x(1)];
elseif P==35
    %'EASY-GEO_PROB'
    J = [-2*x(1) -2*x(2) -2*x(3)];
elseif P==36
    %'EASY-PSS'
    t = Prob.LS.t;
    if(x(3)==(x(1)+x(2)))
        J = [0 0 0 1 0];
    else
        J = zeros(5,5);
        for i=1:5
            J(i,:) =[-x(4)/(x(3)-x(1)-x(2))^2*((x(3)-x(1))*exp(-(x(3)-x(1))*t(i))-x(2)*exp(-x(2)*t(i)))-x(4)/(x(3)-x(1)-x(2))*...
                (-exp(-(x(3)-x(1))*t(i))+(x(3)-x(1))*t(i)*exp(-(x(3)-x(1))*t(i))) ...
                -x(4)/(x(3)-x(1)-x(2))^2*((x(3)-x(1))*exp(-(x(3)-x(1))*t(i))-x(2)*exp(-x(2)*t(i)))-x(4)/(x(3)-x(1)-...
                x(2))*(-exp(-x(2)*t(i))+x(2)*t(i)*exp(-x(2)*t(i))) ...
                x(4)/(x(3)-x(1)-x(2))^2*((x(3)-x(1))*exp(-(x(3)-x(1))*t(i))-x(2)*exp(-x(2)*t(i)))-x(4)/(x(3)-x(1)-...
                x(2))*(exp(-(x(3)-x(1))*t(i))-(x(3)-x(1))*t(i)*exp(-(x(3)-x(1))*t(i))) ...
                1-1/(x(3)-x(1)-x(2))*((x(3)-x(1))*exp(-(x(3)-x(1))*t(i))-x(2)*exp(-x(2)*t(i))) 1];
        end
    end
elseif P==37
    %'EASY-4BAR_LNK'
    t = Prob.LS.t;
    J1 = zeros(24,4);
    J2 = zeros(24,4);
    p1 = zeros(24,1);
    p2 = zeros(24,1);
    a = zeros(24,1);
    d = zeros(24,1);
    for i=1:24
        if((x(1)==0 & x(4)==0) | x(3)==0)%Avoid divide by zero
            p1(i) = 1e9;
            p2(i) = 1e9;
            a(i) = 1e9;
            d(i) = 1e9;
        else
            p1(i) = x(1)*cos(1/12*pi*t(i));
            p2(i) = x(1)*cos(1/12*pi*t(i));
            a(i) = (1/12)*p1(i)*t(i);
            d(i) = x(1)^2+1+2*p1(i);
        end
        J1(i,:) = [-2*(1/2*(2*x(1)+2*cos(a(i)))/x(3)/d(i)^(1/2)-1/4*(x(3)^2+d(i)-x(2)^2)/x(3)/d(i)^(3/2)*(2*x(1)+...
            2*cos(a(i))))/(4-(x(3)^2+d(i)-x(2)^2)^2/x(3)^2/d(i))^(1/2)+...
            (sin(a(i))/(x(4)+p1(i))-p2(i)/(x(4)+p1(i))^2*cos(a(i)))/(1+x(1)^2*sin(a(i))^2/(x(4)+p1(i))^2) ...
            2*x(2)/x(3)/d(i)^(1/2)/(4-(x(3)^2+d(i)-x(2)^2)^2/x(3)^2/d(i))^(1/2) ...
            -2*(1/(sqrt(d(i)))-1/2*(x(3)^2+d(i)-x(2)^2)/x(3)^2/d(i)^(1/2))/(4-(x(3)^2+d(i)-x(2)^2)^2/x(3)^2/d(i))^(1/2) ...
            -p2(i)/(x(4)+p1(i))^2/(1+x(1)^2*sin(a(i))^2/(x(4)+p1(i))^2)];
    end
    J = J1-J2;
elseif P==38
    %'EASY-POL_APP'
    t = Prob.LS.t;
    J = zeros(19,14);
    for i=1:19
        if(t(i) <= 90)
            J(i,:) = [1 t(i) t(i)^2 t(i)^3 t(i)^4 t(i)^5 t(i)^6 0 0 0 0 0 0 0];
        else
            J(i,:) = [0 0 0 0 0 0 0 1 t(i)-90 (t(i)-90)^2 (t(i)-90)^3 (t(i)-90)^4 (t(i)-90)^5 (t(i)-90)^6];
        end
    end
elseif P==39
    %'EASY-LIN_CMP1'
    J = zeros(11,7);
    t = Prob.LS.t;
    p1 = zeros(11,1);
    p2 = zeros(11,1);
    for i=1:11
        p1(i) = exp(x(2)*(t(i)-x(7)));
        p2(i) = exp(x(3)*(t(i)-x(7)));
        if (t(i) >= x(7))
            if (x(2) ~= x(3))
                a2=-(x(5)+x(2))/(x(2)-x(3));
                a1=-1-a2;
                J(i,:) = [a1*p1(i)+a2*p1(i)+1 ...
                    x(1)*(((1/(x(2)-x(3)))+a2)*p1(i)+(a1)*(t(i)-x(7))*p1(i)-p2(i)/(x(2)-x(3))+((x(5)+x(2))*p2(i)/(x(2)-x(3))^2)) ...
                    x(1)*((-(a2)*p1(i)/(x(2)-x(3)))-(-(a2)*p2(i)/(x(2)-x(3)))-(-(a2)*p2(i))) ...
                    0 ...
                    x(1)*((p1(i)/(x(2)-x(3)))-(p2(i)/(x(2)-x(3)))) ...
                    0 ...
                    x(1)*(-a1*x(2)*p1(i)-a2*x(3)*p2(i))];
            else
                J(i,:)=[-p1(i)+1 -x(1)*(t(i)-x(7))*p1(i) 0 0 0 0 x(1)*x(2)*p1(i)];
            end
        else
            J = zeros(11,7);
        end
    end
elseif P==40
    %'EASY-DNS'
    J = zeros(30,3);
    t = Prob.LS.t;
    for i=1:30
        J(i,:) = [x(2)*(exp(-x(3)*t(i))-exp(-x(2)*t(i)))/(x(2)-x(3)) ...
            x(1)*(exp(-x(3)*t(i))-exp(-x(2)*t(i)))/(x(2)-x(3))+x(1)*x(2)*t(i)*exp(-x(2)*t(i))/(x(2)-x(3))-...
            x(1)*x(2)*(exp(-x(3)*t(i))-exp(-x(2)*t(i)))/(x(2)-x(3))^2 ...
            -x(1)*x(2)*t(i)*exp(-x(3)*t(i))/(x(2)-x(3))+x(1)*x(2)*(exp(-x(3)*t(i))-exp(-x(2)*t(i)))/(x(2)-x(3))^2];
    end
elseif P==41
    %'EASY-RTD'
    J = zeros(26,2);
    t = Prob.LS.t;
    for i=1:26
        J(i,:) = [1/(x(2)-x(1))^2*(exp(-t(i)/x(2))-exp(-t(i)/x(1)))-1/(x(2)-x(1))*t(i)/x(1)^2*exp(-t(i)/x(1)) ...
            -1/(x(2)-x(1))^2*(exp(-t(i)/x(2))-exp(-t(i)/x(1)))+1/(x(2)-x(1))*t(i)/x(2)^2*exp(-t(i)/x(2))];
    end
elseif P==42
    %'EASY-RAT_APP'
    J = zeros(11,4);
    t = Prob.LS.t;
    for i=1:11
        J(i,:) = [(t(i)^2+x(2)*t(i))/(t(i)^2+x(3)*t(i)+x(4))  x(1)*t(i)/(t(i)^2+x(3)*t(i)+x(4)) ...
            -x(1)*(t(i)^2+x(2)*t(i))/(t(i)^2+x(3)*t(i)+x(4))^2*t(i)  -x(1)*(t(i)^2+x(2)*t(i))/(t(i)^2+x(3)*t(i)+x(4))^2];
    end
elseif P==43
    %'EASY-MAX_PAT1'
    J = zeros(38,2);
    t = Prob.LS.t;
    for i = 1:38
        if (x(1) ~= x(2) & x(2) ~= 0)
            J(i,:) = [1/(x(2)-x(1))^2*(exp(-t(i)/x(2))-exp(-t(i)/x(1)))-1/(x(2)-x(1))*t(i)/x(1)^2*exp(-t(i)/x(1)) ...
                -1/(x(2)-x(1))^2*(exp(-t(i)/x(2))-exp(-t(i)/x(1)))+1/(x(2)-x(1))*t(i)/x(2)^2*exp(-t(i)/x(2))];
        else
            J = zeros(38,2);
        end
    end
elseif P==44
    %'EASY-MAX_PAT2'
    J = zeros(34,3);
    t = Prob.LS.t;
    for i = 1:34
        if (x(2) ~= x(3)  & x(3) ~= 0 & x(2) ~= 0)
            J(i,:) = [1/(x(3)-x(2))*(exp(-t(i)/x(3))-exp(-t(i)/x(2)))-1/x(3)*exp(-t(i)/x(3)) ...
                x(1)/(x(3)-x(2))^2*(exp(-t(i)/x(3))-exp(-t(i)/x(2)))-x(1)/(x(3)-x(2))*t(i)/x(2)^2*exp(-t(i)/x(2)) ...
                -x(1)/(x(3)-x(2))^2*(exp(-t(i)/x(3))-exp(-t(i)/x(2)))+x(1)/(x(3)-x(2))*t(i)/x(3)^2*exp(-t(i)/x(3))-...
                (1-x(1))/x(3)^2*exp(-t(i)/x(3))+(1-x(1))/x(3)^3*t(i)*exp(-t(i)/x(3))];
        else
            J = zeros(34,3);
        end
    end
elseif P==45
    %'EASY-TREND'
    a = 971;
    t = Prob.LS.t;
    J = zeros(500,6);
    for i=1:500
        J(500,:)=[1  -(a-t(i))^x(5)  -(a-t(i))^x(5)*cos(x(6)*log10(a-t(i)))  (a-t(i))^x(5)*sin(x(6)*log10(a-t(i)))...
            -x(2)*(a-t(i))^x(5)*log10(a-t(i))-x(3)*(a-t(i))^x(5)*log10(a-t(i))*cos(x(6)*log10(a-t(i)))+...
            x(4)*(a-t(i))^x(5)*log10(a-t(i))*sin(x(6)*log10(a-t(i))) ...
            x(3)*(a-t(i))^x(5)*sin(x(6)*log10(a-t(i)))*log10(a-t(i))+x(4)*(a-t(i))^x(5)*cos(x(6)*log10(a-t(i)))*log10(a-t(i))];
    end
end

% MODIFICATION LOG:
%
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981027  hkh  Jacobian for P==17 corrected.
% 981210  hkh  Add constrained Walsh problem. Enables convergence with worse x_0
% 020110  hkh  Add safe guard on Walsh problem (P==28) because of MINOS
% 040517  med  Added 18 problems
% 080603  med  Switched to clsAssign, cleaned