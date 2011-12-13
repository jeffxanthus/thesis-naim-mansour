% function dc=cls_dc(x, Prob)
%
% The cls problems are Constrained Nonlinear Least Squares Problems
%
% cls computes the gradient to the nonlinear constraints c in the point x for
% the test problem P (Prob.P)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function dc=cls_dc(x, Prob)

P=Prob.P;
dc=[];

if P==1
    % 'DIST'
elseif P==2
    % 'BAZA'
elseif P==3
    % 'GMW'
elseif P==4
    % 'PAW'
elseif P==5
    % 'DAS 1'
elseif P==6
    % 'DAS 2'
elseif P==7
    % 'DAS 3'
elseif P==8
    % 'DAS 4'
elseif P==9
    % 'Bob8'
elseif P==10
    % 'Bob9'
elseif P==11
    % 'TP001'
elseif P==12
    % 'TP002'
elseif P==13
    % 'TP028'
elseif P==14
    % 'TP032'
elseif P==15
    % 'TP048'
elseif P==16
    % 'TP049'
elseif P==17
    % 'TP050'
elseif P==18
    % 'TP051'
elseif P==19
    % 'TP052'
elseif P==20
    % 'TP053'
elseif P==21
    % 'TP224'
elseif P==22
    % 'TP231'
elseif P==23
    % 'TP269'
elseif P==24
    % 'TP354'
elseif (P==25)|(P==26)
    % 'WrHo1 and WrHo2';
elseif P==27
    % 'RELN'
elseif P==28
    % Constrained Walsh
elseif P==29
    %'EASY-TP14'
    dc=[-0.5*x(1) -2*x(2)];
elseif P==30
    %'EASY-TP6'
    dc = [-20*x(1) 10];
elseif P==31
    %'EASY-TP43'
    dc = [-2*x(1)-1  -2*x(2)+1  -2*x(3)-1 -2*x(4)+1;
        -2*x(1)+1  -4*x(2)    -2*x(3) 1;
        -4*x(1)-2  -2*x(2)+1  -2*x(3) 1];
elseif P==32
    %'EASY-TP57'
    dc = [-x(2) 0.49];
elseif P==33
    %'EASY-TP327'
    dc = [-x(2) 0.49];
elseif P==34
    % EASY-TP355'
    dc1 = 30*x(2)^2*x(4)^2*x(3)+60*x(2)*x(4)*x(1)-12*x(2)*x(4)*x(3)-126*x(2)^2*x(4)^2*x(1)-6*x(1)...
        +60*x(4)-60*x(2)-30*x(1)*x(4)^2+6*x(3)*x(4)^2+300*x(2)^2*x(4);
    dc2 = 600*x(2)*x(4)*x(1)-6*x(2)*x(4)^2*x(3)^2-12*x(1)*x(3)*x(4)-120*x(2)*x(4)*x(3)-126*x(2)*x(4)^2*x(1)^2 ...
        -60*x(1)+6*x(4)-600*x(2)-30*x(2)*x(4)^2+30*x(1)^2*x(4)+60*x(2)*x(4)^2*x(3)*x(1);
    dc3 = -6*x(2)^2*x(4)^2*x(3)-12*x(2)*x(4)*x(1)+30*x(2)^2*x(4)^2*x(1)+1.500000000*x(3)+6*x(1)*x(4)^2-60*x(2)^2*x(4);
    dc4 = -6*x(2)^2*x(4)*x(3)^2+12*x(1)*x(3)*x(4)-126*x(2)^2*x(4)*x(1)^2-12*x(1)*x(2)*x(3)+60*x(1)-6*x(4)+6*x(2)...
        -30*x(2)^2*x(4)-30*x(1)^2*x(4)+60*x(2)^2*x(4)*x(3)*x(1)+30*x(1)^2*x(2)-60*x(2)^2*x(3)+300*x(1)*x(2)^2;
    dc = [dc1 dc2 dc3 dc4];
elseif P==35
    %'EASY-GEO_PROB'
    dc = [8*x(1)-4+0.1*x(2)  4*x(2)-0.8+0.1*x(1)+0.2*x(3) 2*x(3)+0.2*x(2);
        4*x(1) 2*x(2) -4*x(3)];
elseif P==36
    %'EASY-PSS'
elseif P==37
    %'EASY-4BAR_LNK'
elseif P==38
    %'EASY-POL_APP'
    t = Prob.LS.t;
    dc = zeros(19,14);
    for i = 1:19
        dc(i,:)=[1 t(i) t(i)^2 t(i)^3 t(i)^4 t(i)^5 t(i)^6 -1 0 0 0 0 0 0];
    end
elseif P==39
    % EASY-LIN_CMP1'
    dc = [0 x(3) x(2) 0 -x(6) -x(5) 0];
elseif P==40
    %'EASY-DNS'
elseif P==41
    %'EASY-RTD'
elseif P==42
    %'EASY-RAT_APP'
    TE = 4.0;
    T1 = 0.0625;
    dc = [T1*(T1+x(2))/(T1^2+x(3)*T1+x(4)) T1*x(1)/(T1^2+x(3)*T1+x(4))   -x(1)*T1^2*(T1+x(2))/(T1^2+x(3)*T1+x(4))^2 -x(1)*T1*(T1+x(2))/(T1^2+x(3)*T1+x(4))^2;
        TE*(TE+x(2))/(TE^2+x(3)*TE+x(4))   TE*x(1)/(TE^2+x(3)*TE+x(4))   -x(1)*TE^2*(TE+x(2))/(TE^2+x(3)*TE+x(4))^2  -x(1)*TE*(TE+x(2))/(TE^2+x(3)*TE+x(4))^2];
elseif P==43
    %'EASY-MAX_PAT1'
elseif P==44
    %'EASY-MAX_PAT2'
elseif P==45
    %'EASY-TREND'
    dc = [0 2*x(2) -2*x(3)*(x(6)^2+x(5)^2)/x(5)^2 -2*x(4)*(x(6)^2+x(5)^2)/x(5)^2 -2*(x(3)^2+x(4)^2)/x(5)+2*(x(3)^2+x(4)^2)*(x(6)^2+x(5)^2)/x(5)^3 -2*(x(3)^2+x(4)^2)*x(6)/x(5)^2];
end

% MODIFICATION LOG:
%
% 040517  med  Added 18 problems
% 051128  bjo  Constraint Jacobian for P==42 corrected
% 080603  med  Switched to clsAssign, cleaned