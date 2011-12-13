% goals_r.m
%
% function r = goals_r(x, Prob)
%
% Computes residuals for
% Multi Criterium Unconstrained & Constrained Nonlinear Problem
% in the point x for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function r = goals_r(x, Prob)

P = Prob.P;
y = Prob.LS.y;
m = size(y,1);

if P==1
    % EASY-FIT 'TP269'
    r= [x(2)-x(1); x(2)+x(3)-2; x(4)-1; x(5)-1];
elseif P==2
    % EASY-FIT 'TP13'
    r=[2-x(1);x(2)];
elseif P==3
    % EASY-FIT 'TP46'
    r=[x(1)-x(2); x(3)-1;(x(4)-1)^2;(x(5)-1)^3];
elseif P==4
    % EASY-FIT 'TP48'
    r=[x(1)-1; x(2)-x(3); x(4)-x(5)];
elseif P==5
    % EASY-FIT 'TP354'
    r=[x(1)+10*x(2);sqrt(5)*(x(3)-x(4)); (x(2)-2*x(3))^2; sqrt(10)*(x(1)-x(4))^2];
elseif P==6
    % EASY-FIT 'TP355'
    r=[11-x(1)*x(4)-x(2)*x(4)+x(3)*x(4); x(1)+10*x(2)-x(3)+x(4)+x(2)*x(4)*(x(3)-x(1))];
elseif P==7
    % EASY-FIT 'TP372'
    r=[x(4);x(5);x(6);x(7);x(8);x(9)];
elseif P==8
    % EASY-FIT 'TP373'
    r=[x(4);x(5);x(6);x(7);x(8);x(9)];
elseif P==9
    %'EASY-OPT_KINX'
    t = Prob.LS.t;
    k12=0.1;
    k23=0.055;
    V2=7.0;
    y2_x4 = x(3)/k23*(1 - exp(-k23*x(4)));
    y2_x5=y2_x4*exp(-k23*(x(5)-x(4)));
    a1 = y2_x5 - k12*x(1)/(k23 - k12);
    b1 = k12*x(1)/(k23 - k12);
    y2_x6 = a1*exp(-k23*(x(6)-x(5)))+b1*exp(-k12*(x(6)-x(5)));
    a2 = y2_x6 - k12*x(2)/(k23 - k12);
    b2 = k12*x(2)/(k23 - k12);
    r=zeros(100,1);
    for i=1:100
        if(t(i) <= x(4))
            r(i) = x(3)/k23*(1-exp(-k23*t(i)))/V2;
        elseif(t(i) <=x(5))
            r(i) = y2_x4*exp(-k23*(t(i)-x(4)))/V2;
        elseif(t(i)<=x(6))
            r(i) = (a1*exp(-k23*(t(i)-x(5))) + b1*exp(-k12*(t(i)-x(5))))/V2;
        else
            r(i) = (a2*exp(-k23*(t(i)-x(6))) + b2*exp(-k12*(t(i)-x(6))))/V2;
        end
    end
end

if Prob.LS.yUse & m==length(r),  r=r-y; end

% MODIFICATION LOG:
%
% 040517  med  Created
% 050508  med  Print of r for problem 1 removed
% 080603  med  Switched to clsAssign, cleaned