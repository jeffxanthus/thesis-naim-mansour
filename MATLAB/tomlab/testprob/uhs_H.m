% function H = uhs_H(x,Prob)
%
% Hessian matrix for test functions for Unconstrained optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function H = uhs_H(x,Prob)

x=x(:);
P=Prob.P;

if P == 1 % HS 1
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)
        -400*x(1)             200   ];
elseif P == 2 % HS 2
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)
        -400*x(1)             200   ];
elseif P == 3 % HS 3
    H = 1/50000*[1 -1;-1  1];
elseif P == 4 % HS 4
    H = [2*x(1)+2 0;0 0];
elseif P == 5 % HS 5
    H = [-sin(x(1)+x(2))+2 -sin(x(1)+x(2))-2
        -sin(x(1)+x(2))-2 -sin(x(1)+x(2))+2];
elseif P == 6 % HS 38
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)  0                      0
        -400*x(1)               220.2      0                      19.8
        0                       0          1080*x(3)^2-360*x(4)+2 -360*x(3)
        0                       19.8       -360*x(3)              200.2];
elseif P == 7 % HS 45
    H = -1/120 * [...
        0,x(3)*x(4)*x(5),x(2)*x(4)*x(5),x(2)*x(3)*x(5),x(2)*x(3)*x(4) ; ...
        0,0,x(1)*x(4)*x(5),x(1)*x(3)*x(5),x(1)*x(3)*x(4) ; ...
        0,0,0,x(1)*x(2)*x(5)*x(5),x(1)*x(2)*x(4) ; ...
        0,0,0,0,x(1)*x(2)*x(3) ; ...
        0,0,0,0,0 ];
    H = 0.5*(H+H');
elseif P == 8 % HS 110
    n = 10;
    ep = .2*(x(1)*x(2)*x(3)*x(4)*x(5)*x(6)*x(7)*x(8)*x(9)*x(10))^0.2;
    z = zeros(n,1);
    for i = 1:n
        z(i)=2*1/((x(i)-2)^2)-2*log(x(i)-2)/((x(i)-2)^2) ...
            +2*1/((10-x(i))^2)-2*log(10-x(i))/((10-x(i))^2);
    end
    H = diag(z)-(diag(-5*ones(n,1))+ones(n,n)).*0.2*ep./(x*x');
elseif P == 9 % HS 206
    H = [12*x(1)^2-4*x(2)+200  -4*x(1)
        -4*x(1)          2     ];
elseif P == 10 % HS 207
    H = [12*x(1)^2-4*x(2)+2    -4*x(1)
        -4*x(1)          2     ];
elseif P == 11 % HS 208
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)
        -400*x(1)             200   ];
elseif P == 12 % HS 209
    H = [120000*x(1)^2-40000*x(2)+2  -40000*x(1)
        -40000*x(1)             20000   ];
elseif P == 13 % HS 210
    H = [12000000*x(1)^2-4000000*x(2)+2  -4000000*x(1)
        -4000000*x(1)             2000000   ];
elseif P == 14 % HS 211
    H = [1800*x(1)^4-1200*(x(2)-x(1)^3)*x(1)+2 -600*x(1)^2
        -600*x(1)^2                             200      ];
elseif P == 15 % HS 212
    H=zeros(2,2);
    H(1,1)=(32+2*(3+(x(1)-2)^2+x(2)^2+(x(1) ...
        -x(2))*(2*x(1)-4))^2+2*(4*x(1)+  ...
        4*x(2)+(x(1)-x(2))*((x(1)-2)^2+  ...
        x(2)^2-1))*(6*x(1)-8-2*x(2)));
    H(1,2)=32+2*(3+(x(1)-2)^2+x(2)^2+(x(1)-x(2))*(2*   ...
        x(1)-4))*(5-(x(1)-2)^2-x(2)^2+2*(x(1)-x(2)) ...
        *x(2))+2*(4*x(1)+4*x(2)+(x(1)-x(2))*((x(1)  ...
        -2)^2+x(2)^2-1))*(-2*x(1)+4+2*x(2));
    H(2,1)=H(1,2);
    H(2,2)=32+2*(5-(x(1)-2)^2-x(2)^2+2*(x(1)-x(2))*x(2))^2+2*(4*x(1)+ ...
        4*x(2)+(x(1)-x(2))*((x(1)-2)^2+x(2)^2-1))*(-6*x(2)+2*x(1));
elseif P == 16 % HS 213
    H=zeros(2,2);
    H(1,1)=12*(10*(x(1)-x(2))^2+     ...
        (-1+x(1))^2)^2*(22*x(1)-  ...
        20*x(2)-2)^2+88*(10*(x(1) ...
        -x(2))^2+(-1+x(1))^2)^3;
    H(1,2)=12*(10*(x(1)-x(2))^2+(-1+     ...
        x(1))^2)^2*(-20*x(1)+20*x(2)) ...
        *(22*x(1)-20*x(2)-2)-80*(10*  ...
        (x(1)-x(2))^2+(-1+x(1))^2)^3;
    H(2,1)=H(1,2);
    H(2,2)=12*(10*(x(1)-x(2))^2+(-1+x(1))          ...
        ^2)^2*(-20*x(1)+20*x(2))^2+80*(10*(x(1) ...
        -x(2))^2+(-1+x(1))^2)^3;
elseif P == 17 % HS 214
    H=zeros(2,2);
    H(1,1)=-3/16*(22*x(1)-20*x(2)-2)   ...
        ^2/((10*(x(1)-x(2))^2+(-1+  ...
        x(1))^2)^(7/4))+11/2*1/((10 ...
        *(x(1)-x(2))^2+(-1+x(1))^2)^(3/4));
    H(1,2)=-3/16*(-20*x(1)+20*x(2))     ...
        *(22*x(1)-20*x(2)-2)/((10*   ...
        (x(1)-x(2))^2+(-1+x(1))^2)   ...
        ^(7/4))-5*1/((10*(x(1)-x(2)) ...
        ^2+(-1+x(1))^2)^(3/4));
    H(2,1)=H(1,2);
    H(2,2)=-3/16*(-20*x(1)+20*x(2))  ...
        ^2/((10*(x(1)-x(2))^2+(-1 ...
        +x(1))^2)^(7/4))+5*1/((10 ...
        *(x(1)-x(2))^2+(-1+x(1))^2)^(3/4));
elseif P == 18 % HS 229
    H = [ 1200*x(1)^2-400*x(2)+2  -400*x(1)
        -400*x(1)               200        ];
elseif P == 19 % HS 247
    if x(1) > 0
        u = 0.5*pi*atan(x(2)/x(1));
    elseif x(1) < 0
        u = 0.5+0.5*pi*atan(x(2)/x(1));
    else
        u = inf;
    end
    H=zeros(3,3);
    H(1,1)=5000*1/(pi^2*x(2)^2*(1+x(1)^2/(x(2)^2))^2) ...
        +2000*(x(3)-10*u)*x(1)/(pi*x(2)^3*(1+x(1)^2 ...
        /(x(2)^2))^2)+200*x(1)^2/(x(1)^2+x(2)^2)-200 ...
        *(sqrt(x(1)^2+x(2)^2)-1)*x(1)^2/((x(1)^2+x(2) ...
        ^2)^(3/2))+200*(sqrt(x(1)^2+x(2)^2)-1)/(sqrt(x ...
        (1)^2+x(2)^2));
    H(1,2)=-5000*x(1)/(pi^2*x(2)^3*(1+x(1)^2/(x(2)^2))^2) ...
        +1000*(x(3)-10*u)/(pi*x(2)^2*(1+x(1)^2/(x(2) ...
        ^2)))-2000*(x(3)-10*u)*x(1)^2/(pi*x(2)^4*(1+x(1)^2/ ...
        (x(2)^2))^2)+200*x(2)*x(1)/(x(1)^2+x(2)^2)-200*(sqrt ...
        (x(1)^2+x(2)^2)-1)*x(1)*x(2)/((x(1)^2+x(2)^2)^(3/2));
    H(1,3)=-1000*1/(pi*x(2)*(1+x(1)^2/(x(2)^2)));
    H(2,1)=H(1,2);
    H(2,2)=5000*x(1)^2/(pi^2*x(2)^4*(1+x(1)^2/(x(2)^2) ...
        )^2)-2000*(x(3)-10*u)*x(1)/(pi*x(2)^3*(1+x(1) ...
        ^2/(x(2)^2)))+2000*(x(3)-10*u)*x(1)^3/(pi*x(2 ...
        )^5*(1+x(1)^2/(x(2)^2))^2)+200*x(2)^2/(x(1)^2 ...
        +x(2)^2)-200*(sqrt(x(1)^2+x(2)^2)-1)*x(2)^2/( ...
        (x(1)^2+x(2)^2)^(3/2))+200*(sqrt(x(1) ...
        ^2+x(2)^2)-1)/(sqrt(x(1)^2+x(2)^2));
    H(2,3)=1000*x(1)/(pi*x(2)^2*(1+x(1)^2/(x(2)^2)));
    H(3,1)=H(1,3);
    H(3,2)=H(2,3);
    H(3,3)=202;
elseif P == 20 % HS 256
    H = [2+120*(x(1)-x(4))^2  20  0  -120*(x(1)-x(4))^2
        20  200+12*(x(2)-2*x(3))^2  -24*(x(2)-2*x(3))^2  0
        0  -24*(x(2)-2*x(3))^2  10+48*(x(2)-2*x(3))^2  -10
        -120*(x(1)-x(4))^2  0  -10  10+120*(x(1)-x(4))^2];
elseif P == 21 % HS 257
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)  0                       19.8
        -400*x(1)               220.2      0                       0
        0                       0          1080*x(3)^2-360*x(4)+2  -360*x(3)
        19.8                    0          -360*x(3)               200.2];
elseif P == 22 % HS 258
    H = [1200*x(1)^2-400*x(2)+2  -400*x(1)  0                       0
        -400*x(1)               220.2      0                       19.8
        0                       0          1080*x(3)^2-360*x(4)+2  -360*x(3)
        0                       19.8       -360*x(3)               200.2];
elseif P == 23 % HS 259
    H = [800*x(1)^2-400*(-x(1)^2+x(2))+2  -400*x(1)  0              0
        -400*x(1)               220.2      0                       19.8
        0            0          -360*(-x(3)^2+x(4))-360*(-2*x(3))*x(3)+6*(1-x(3))  -360*x(3)
        0            19.8       -360*x(3)                          182];
elseif P == 24 % HS 260
    H = [1200*x(1)^92-400*x(2)+2  -400*x(1)  0  0
        -400*x(1)  220.2  0  19.8
        0  0  1080*x(3)^2-360*x(4)+2  -360*x(3)
        0 19.8 -360*x(3) 200.2];
elseif P == 25 % HS 273
    H = [300+20*(30*x(1)-30)^2+9000*(x(1)-1)^2+8400*(x(2)-1) ...
        ^2+7800*(x(3)-1)^2+7200*(x(4)-1)^2+6600*(x(5)-1)^2+6000 ...
        *(x(6)-1)^2, 20*(30*x(1)-30)*(28*x(2)-28), 20*(30 ...
        *x(1)-30)*(26*x(3)-26), 20*(30*x(1)-30)*(24*x(4)-24), 20 ...
        *(30*x(1)-30)*(22*x(5)-22),  20*(30*x(1)-30)*(20*x(6)-20);
        20*(30*x(1)-30)*(28*x(2)-28), 280+20*(28*x(2)-28)^2+ ...
        8400*(x(1)-1)^2+7840*(x(2)-1)^2+7280*(x(3)-1)^2+6720 ...
        *(x(4)-1)^2+6160*(x(5)-1)^2+5600*(x(6)-1)^2, 20*(28* ...
        x(2)-28)*(26*x(3)-26), 20*(28*x(2)-28)*(24*x(4)-24), ...
        20*(28*x(2)-28)*(22*x(5)-22), 20*(28*x(2)-28)*(20*x(6)-20)];...
        h3=[
        20*(30*x(1)-30)*(26*x(3)-26), 20*(28*x(2)-28)*(26*x(3)-26), ...
        260+20*(26*x(3)-26)^2+7800*(x(1)-1)^2+7280*(x(2)-1)^2+ ...
        6760*(x(3)-1)^2+6240*(x(4)-1)^2+5720*(x(5)-1)^2+ ...
        5200*(x(6)-1)^2, 20*(26*x(3)-26)*(24*x(4)-24),  20*(26* ...
        x(3)-26)*(22*x(5)-22), 20*(26*x(3)-26)*(20*x(6)-20)];
    h4=[20*(30*x(1)-30)*(24*x(4)-24),...
        20*(28*x(2)-28)*(24*x(4)-24),...
        20*(26*x(3)-26)*(24*x(4)-24),...
        240+20*(24*x(4)-24)^2+7200*(x(1)-1)^2+ ...
        6720*(x(2)-1)^2+6240*(x(3)-1)^2+ ...
        5760*(x(4)-1)^2+5280*(x(5)-1)^2+4800*(x(6)-1)^2, ...
        20*(24*x(4)-24)*(22*x(5)-22), ...
        20*(24*x(4)-24)*(20*x(6)-20)];
    h5=[20*(30*x(1)-30)*(22*x(5)-22), 20*(28*x(2)-28)*(22*x(5)-22),...
        20*(26*x(3)-26)*(22*x(5)-22), 20*(24*x(4)-24)*...
        (22*x(5)-22), 220+20*(22*x(5)-22)^2+6600*(x(1)-1)^2+...
        +6160*(x(2)-1)^2+5720*(x(3)-1)^2+5280*(x(4)-1)^2+...
        4840*(x(5)-1)^2+4400*(x(6)-1)^2, 20*(22*x(5)-22)*(20*x(6)-20)];
    h6=[20*(30*x(1)-30)*(20*x(6)-20), 20*(28*x(2)-28)*(20*x( ...
        6)-20), 20*(26*x(3)-26)*(20*x(6)-20), 20*(24*x(4)-24 ...
        )*(20*x(6)-20), 20*(22*x(5)-22)*(20*x(6)-20), 200+20 ...
        *(20*x(6)-20)^2+6000*(x(1)-1)^2+5600*(x(2)-1)^2+5200 ...
        *(x(3)-1)^2+4800*(x(4)-1)^2+4400*(x(5)-1)^2+4000*(x(6)-1)^2];
    H = [H;h3;h4;h5;h6];
end

% MODIFICATION LOG:
%
% 990601  mbk  First version.
% 030515  ango Fixed error in Problem 7
% 050215  med  Updated problem 23
% 080603  med  Switched to conAssign, cleaned