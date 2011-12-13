% mco_J.m
%
% function J = mco_J(x, Prob)
%
% Computes the Jacobian to the Unconstrained & Constrained Nonlinear Multi Criterium Problem
% in the point x for the test problem P (Prob.P).
%
% J(i,j) is dr_i/d_x_j

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function J = mco_J(x, Prob)

P=Prob.P;

if P==1
    % 'MC0-TP 1'
    p1=(x(2)-0.2)/0.004;
    p2=(x(2)-0.6)/0.4;
    g=-2.0-exp(-(p1)^2)-0.8*exp(-(p2)^2);
    J=[1 0;
        -g/(x(1))^2 (-(-125000*x(2)+25000)*exp(-(p1)^2)-0.8*(-12.5*x(2)+7.5)*exp(-(p2)^2))/x(1)];
elseif P==2
    % 'MCO-TP 2'
    if x(1) <= 1
        J = [-1 ;2*(x(1)-5) ];
    elseif x(1) >1 & x(1)<=3
        J = [ 1 ;2*(x(1)-5) ];
    elseif x(1) >3 & x(1)<=4
        J = [-1 ;2*(x(1)-5) ];
    else % if x(1)>4
        J = [ 1 ;2*(x(1)-5) ];
    end
elseif P==3
    %'MCO-TP 3'
    g = 11+x(2)^2-10*cos(2*pi*x(2));
    if(x(1)==0)
        dr2x1=-1e9;% Obs!  OPTIONAL BOUND LARGE TO AVOID DIVIDE BY ZERO
    else
        dr2x1=-(1/2)*sqrt(g/x(1));
    end
    if(x(1)<=g)
        J = [1 0; dr2x1 2*x(2)+20*pi*sin(2*pi*x(2))*(1-0.5*sqrt(x(1)/g))];
    else
        J =  [1 0;0 0];
    end
elseif P==4
    %'MCO-TP 4'
    p1=(x(1)-1)^2+(x(2)-0.5)^2;
    p2=(x(1)+1)^2+(x(2)+0.5)^2;
    p3=(x(1)-1)^2+(x(2)+0.5)^2;
    J = [2*(x(1)-1)*exp(-p1) 2*(x(2)-0.5)*exp(-p1);
        2*(x(1)+1)*exp(-p2) 2*(x(2)+0.5)*exp(-p2);
        2*(x(1)-1)*exp(-p3) 2*(x(2)+0.5)*exp(-p3)];
elseif P==5
    %'MCO-TP 5'
    sum = 0;
    for i = 2:10
        sum = sum + x(i)^2-10*cos(4*pi*x(i));
    end
    g = 91 + sum;
    if(x(1)==0)
        dr2x1=-1e9;% Obs!  OPTIONAL BOUND LARGE TO AVOID DIVIDE BY ZERO
    else
        dr2x1=-(1/2)*sqrt(g/x(1));
    end
    J = [1 0 0 0 0 0 0 0 0 0;
        dr2x1 ...
        (2*x(2)+40*sin(4*pi*x(2))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(3)+40*sin(4*pi*x(3))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(4)+40*sin(4*pi*x(4))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(5)+40*sin(4*pi*x(5))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(6)+40*sin(4*pi*x(6))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(7)+40*sin(4*pi*x(7))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(8)+40*sin(4*pi*x(8))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(9)+40*sin(4*pi*x(9))*pi)*(1-(1/2)*sqrt(x(1)/(g))) ...
        (2*x(10)+40*sin(4*pi*x(10))*pi)*(1-(1/2)*sqrt(x(1)/(g)))];
elseif P==6
    %'MCO-TP 6'
    sum = 0;
    prod = 1;
    for i=2:10;
        sum =  sum + (x(i)^2)/4000;
        prod = prod*cos(x(i)/sqrt(i));
        g = 2 + sum - prod;
    end
    if(x(1)==0)
        dr2x1=-1e9;% Obs! OPTIONAL BOUND LARGE TO AVOID DIVIDE BY ZERO
    else
        dr2x1=-(1/2)*sqrt(g/x(1));
    end
    J = [1 0 0 0 0 0 0 0 0 0;
        dr2x1 ...
        (2*(x(2)/4000)+(1/sqrt(2))*sin(x(2)/sqrt(2)))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(3)/4000)+(1/sqrt(3))*sin(x(3)/sqrt(3)))*cos(x(2))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(4)/4000)+(1/sqrt(4))*sin(x(4)/sqrt(4)))*cos(x(2))*cos(x(3))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(5)/4000)+(1/sqrt(5))*sin(x(5)/sqrt(5)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(6))*cos(x(7))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(6)/4000)+(1/sqrt(6))*sin(x(6)/sqrt(6)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(7))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(7)/4000)+(1/sqrt(7))*sin(x(7)/sqrt(7)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(8))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(8)/4000)+(1/sqrt(8))*sin(x(8)/sqrt(8)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(9))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(9)/4000)+(1/sqrt(9))*sin(x(9)/sqrt(9)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(8))*cos(x(10))*(1-(1/2)*sqrt(x(1)/g)) ...
        (2*(x(10)/4000)+(1/sqrt(10))*sin(x(10)/sqrt(10)))*cos(x(2))*cos(x(3))*cos(x(4))*cos(x(5))*cos(x(6))*cos(x(7))*cos(x(6))*cos(x(9))*(1-(1/2)*sqrt(x(1)/g))];
elseif P==7
    J = []; %Not given
elseif P==8
    %'MCO-TP 8'
    J = [2*(x(1)-2) 2*(x(2)-1);
        9 -2*(x(2)-1)];
elseif P==9
    %'MCO-TP 9'
    J = [-50*(x(1)-2) -2*(x(2)-2) -2*(x(3)-1) -2*(x(4)-4) -2*(x(5)-1) 0;
        2*x(1) 2*x(2) 2*x(3) 2*x(4) 2*x(5) 2*x(6)];
end

% MODIFICATION LOG:
%
% 040511  med  Created.
% 050503  hkh  Corrected problem 2, probably still wrong, same with 7
% 080603  med  Switched to clsAssign, cleaned