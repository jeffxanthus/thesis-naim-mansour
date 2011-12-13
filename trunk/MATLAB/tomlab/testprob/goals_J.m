% goals_J.m
%
% function J = goals_J(x, Prob)
%
% Computes the Jacobian for the
% Multi Criterium Unconstrained & Constrained Nonlinear Problem
% in the point x for the test problem P (Prob.P).
%
% J(i,j) is dr_i/d_x_j

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function J = goals_J(x, Prob)

P=Prob.P;

if P==1
    % EASY-FIT 'TP269'
    J = [-1 1 0 0 0;0 1 1 0 0;0 0 0 1 0;0 0 0 0 1];
elseif P==2
    % EASY-FIT 'TP13'
    J = [-1 0; 0 1];
elseif P==3
    % EASY-FIT 'TP46'
    J = [1 -1 0 0 0;0 0 1 0 0;0 0 0 2*(x(4)-1) 0;0 0 0 0 3*(x(5)-1)^2];
elseif P==4
    % EASY-FIT 'TP48'
    J = [1 0 0 0 0;0 1 -1 0 0;0 0 0 1 -1];
elseif P==5
    % EASY-FIT 'TP354'
    J = [1 10 0 0;
        0 0 sqrt(5) -sqrt(5);
        0 2*(x(2)-2*x(3)) -4*(x(2)-2*x(3)) 0;
        2*sqrt(10)*(x(1)-x(4)) 0 0 -2*sqrt(10)*(x(1)-x(4))];
elseif P==6
    J = [-x(4)  -x(4)  x(4) -x(1)-x(2)+x(3);
        1-x(2)*x(4)  10+x(4)*x(3)-x(1)*x(4) -1+x(2)*x(4)  1+x(2)*x(3)-x(2)*x(1)];
elseif P==7
    % EASY-FIT 'TP372'
    J = [0 0 0 1 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 0 1];
elseif P==8
    % EASY-FIT 'TP373'
    J = [0 0 0 1 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 1 0 0 0;
        0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 0 1];
elseif P==9
    %'EASY-OPT_KINX'
    t = Prob.LS.t;
    k12=0.1;
    k23=0.055;
    kr=k12/(k23-k12);
    v2=7.0;
    J =zeros(100,6);
    for i=1:100
        if(t(i) <= x(4))
            J(i,:) = [0 0 2.597402597-2.597402597*exp(-k23*t(i)) 0 0 0];
        elseif(t(i) <= x(5))
            J(i,:) = [0 0 2.597402597*(1-exp(-k23*x(4)))*exp(-k23*t(i)+k23*x(4)) ...
                (1/v2)*x(3)*exp(-k23*x(4))*exp(-k23*t(i)+k23*x(4))+0.1428571428*x(3)*(1-exp(-k23*x(4)))*exp(-k23*t(i)+...
                k23*x(4)) 0 0];
        elseif(t(i)<=x(6))
            J(i,:) = [0.3174603175*exp(-k23*t(i)+k23*x(5))-0.3174603175*exp(-0.1*t(i)+.1*x(5)) 0 2.597402598*(1-...
                exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))*exp(-k23*t(i)+k23*x(5)) ...
                (1/v2)*(0.9999999999*x(3)*exp(-k23*x(4))*exp(-k23*x(5)+k23*x(4))+0.9999999999*x(3)*(1-...
                exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4)))*exp(-k23*t(i)+k23*x(5)) ...
                -(1/v2)*x(3)*(1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))*exp(-k23*t(i)+k23*x(5))+...
                0.7857142860e-2*(18.18181818*x(3)*(1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))+...
                kr*x(1))*exp(-k23*t(i)+k23*x(5))-0.3174603175e-1*x(1)*exp(-0.1*t(i)+0.1*x(5)) 0];
        else
            J(i,:)= [(1/v2)*(kr*exp(-k23*x(6)+k23*x(5))-kr*exp(-0.1*x(6)+0.1*x(5)))*exp(-k23*t(i)+k23*x(6)) ...
                0.3174603175*exp(-k23*t(i)+0.55e-1*x(6))-0.3174603175*exp(-0.1*t(i)+0.1*x(6)) ...
                2.597402598*(1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))*exp(-k23*x(6)+k23*x(5))*exp(-k23*t(i)+k23*x(6)) ...
                0.1428571429*(0.9999999999*x(3)*exp(-k23*x(4))*exp(-k23*x(5)+k23*x(4))+0.9999999999*x(3)*(1-exp(-k23*x(4)))*...
                exp(-k23*x(5)+k23*x(4)))*exp(-k23*x(6)+k23*x(5))*exp(-k23*t(i)+k23*x(6))  (1/v2)*(-0.9999999999*x(3)*...
                (1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))*exp(-k23*x(6)+k23*x(5))+k23*(18.18181818*x(3)*(1-exp(-k23*...
                x(4)))*exp(-k23*x(5)+k23*x(4))+kr*x(1))*exp(-k23*x(6)+k23*x(5))-0.2222222222*x(1)*exp(-0.1*x(6)+0.1*x(5)))*...
                exp(-k23*t(i)+k23*x(6)) (1/v2)*(-k23*(18.18181818*x(3)*(1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))+...
                kr*x(1))*exp(-k23*x(6)+k23*x(5))+0.2222222222*x(1)*exp(-0.1*x(6)+0.1*x(5)))*exp(-k23*t(i)+k23*x(6))+...
                0.7857142860e-2*((18.18181818*x(3)*(1-exp(-k23*x(4)))*exp(-k23*x(5)+k23*x(4))+kr*x(1))*exp(-k23*x(6)+...
                k23*x(5))-kr*x(1)*exp(-0.1*x(6)+0.1*x(5))+kr*x(2))*exp(-k23*t(i)+k23*x(6))-0.3174603175e-1*x(2)*...
                exp(-0.1*t(i)+0.1*x(6))];
        end
    end
end

% MODIFICATION LOG:
%
% 040517  med  Created
% 080603  med  Switched to clsAssign, cleaned