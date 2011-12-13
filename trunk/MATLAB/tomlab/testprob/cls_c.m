% function cx=cls_c(x, Prob)
%
% The cls problems are Constrained Nonlinear Least Squares Problems
%
% cls_c evaluates the nonlinear constraints for test problem P (Prob.P)
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function cx=cls_c(x, Prob)

P=Prob.P;
cx=[];

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
    % 'EASY-TP14'
    cx = -0.25*x(1)^2-x(2)^2;
elseif P==30
    % 'EASY-TP6'
    cx = 10*(x(2)-x(1)^2);
elseif P==31
    % 'EASY-TP43'
    cx = [-x(1)^2-x(2)^2-x(3)^2-x(4)^2-x(1)+x(2)-x(3)+x(4);
        -x(1)^2-2*x(2)^2-x(3)^2-2*x(4)^2+x(1)+x(4);
        -2*x(1)^2-x(2)^2-x(3)^2-2*x(1)+x(2)+x(4)];
elseif P==32
    %'EASY-TP57'
    cx = 0.49*x(2)-x(1)*x(2);
elseif P==33
    %'EASY-TP327'
    cx = 0.49*x(2)-x(1)*x(2);
elseif P==34
    %'EASY-FIT TP355'
    c1 = 11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4);
    c2 = x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)- x(1));
    c3 = 11-4*x(1)*x(4)-4*x(2)*x(4)+x(3)*x(4);
    c4 = 2*x(1)+20*x(2)-0.5*x(3)+2*x(4)+2*x(2)*x(4)*(x(3)-4*x(1));
    cx = (c1)^2+(c2)^2-(c3)^2-(c4)^2;
elseif P==35
    %'EASY-GEO_PROB'
    cx = [4*x(1)^2-4.0*x(1)+2*x(2)^2-0.8*x(2)+x(3)^2+0.1*x(1)*x(2)+0.2*x(2)*x(3);
        2*x(1)^2+x(2)^2-2*x(3)^2];
elseif P==36
    %'EASY-PSS'
elseif P==37
    %'EASY-4BAR_LNK'
elseif P==38
    %'EASY-POL_APP'
    t = Prob.LS.t;
    cx = zeros(19,1);
    for i=1:19
        cx(i) = -x(8)+x(1)+x(2)*t(i)+x(3)*t(i)^2+x(4)*t(i)^3+x(5)*t(i)^4+x(6)*t(i)^5+x(7)*t(i)^6;
    end
elseif P==39
    %'EASY-LIN_CMP1'
    cx = x(2)*x(3)-x(5)*x(6);
elseif P==40
    %'EASY-FIT DNS'
elseif P==41
    %'EASY-FIT RTD'
elseif P==42
    %'EASY-RAT_APP'
    TE = 4.0;
    T1 = 0.0625;
    cx = [x(1)*T1*(T1 + x(2))/(T1^2 + x(3)*T1 + x(4)); x(1)*TE*(TE + x(2))/(TE^2 + x(3)*TE + x(4))];
elseif P==43
    %'EASY-MAX_PAT1'
elseif P==44
    %'EASY-MAX_PAT2'
elseif P==45
    %'EASY-TREND'
    cx = x(2)^2-(x(3)^2+x(4)^2)*(x(6)^2+x(5)^2)/x(5)^2;
end

% MODIFICATION LOG:
%
% 040517  med  Added 18 problems
% 080603  med  Switched to clsAssign, cleaned