% mco_r.m
%
% function r = mco_r(x, Prob)
%
% Computes residuals to Multi Constrained Nonlinear Problem
% in the point x for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function r = mco_r(x, Prob)

P = Prob.P;
y = Prob.LS.y;
m = size(y,1);

if P==1
    % 'MCO-TP 1'
    p1 = (x(2)-0.2)/0.004;
    p2 = (x(2)-0.6)/0.4;
    g  = -2.0-exp(-(p1)^2)-0.8*exp(-(p2)^2);
    r  = [x(1);g/x(1)];
elseif P==2
    % 'MCO-TP 2'
    if x(1) <= 1
        r = [-x(1);(x(1)-5)^2];
    elseif x(1) >1 & x(1)<=3
        r = [ x(1)-2;(x(1)-5)^2];
    elseif x(1) >3 & x(1)<=4
        r = [3-x(1);(x(1)-5)^2];
    else
        r = [x(1)-4;(x(1)-5)^2];
    end
elseif P==3
    % 'MCO-TP 3'
    g = 11+x(2)^2-10*cos(2*pi*x(2));
    h = 1-sqrt(x(1)/g);
    if (x(1) <= g)
        r = [x(1) g*h]';
    else
        r =  [x(1) 0]';
    end
elseif P==4
    % 'MCO-TP 4'
    p1=(x(1)-1)^2+(x(2)-0.5)^2;
    p2=(x(1)+1)^2+(x(2)+0.5)^2;
    p3=(x(1)-1)^2+(x(2)+0.5)^2;
    r = [1-exp(-(p1));1-exp(-(p2));1-exp(-(p3))];
elseif P==5
    % 'MCO-TP 5'
    sum = 0;
    for i = 2:10
        sum = sum + x(i)^2-10*cos(4*pi*x(i));
    end
    g = 91+sum;
    r = [x(1); g*(1-sqrt(x(1)/g))];
elseif P==6
    % 'MCO-TP 6'
    sum = 0;
    prod = 1;
    for i=2:10;
        sum =  sum + (x(i)^2)/4000;
        prod = prod*cos(x(i)/sqrt(i));
        g = 2 + sum - prod;
    end
    r = [x(1); g*(1-sqrt(x(1)/g))];
elseif P==7
    % 'MCO-TP 7'
    sum = 0;
    for i=2:10
        sum = sum + (x(i)/9)^(0.25);
    end
    g = 1 + 9*sum;
    r = [1-(exp(-4*x(1)))*(sin(6*pi*x(1)))^6;g*(1-(x(1)/g)^2)];
elseif P==8
    % 'MCO-TP 8'
    r = [2+(x(1)-2)^2+(x(2)-1)^2;9*x(1)-(x(2)-1)^2];
elseif P==9
    % 'MCO-TP 9'
    r = [-(25*(x(1)-2)^2+(x(2)-2)^2+(x(3)-1)^2+(x(4)-4)^2+(x(5)-1)^2);
        x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2+x(6)^2];
end
if Prob.LS.yUse & m==length(r),  r=r-y; end

% MODIFICATION LOG:
%
% 040511  med  Created.
% 050503  hkh  Corrected problem 2, probably still wrong, same with 7
% 080603  med  Switched to clsAssign, cleaned