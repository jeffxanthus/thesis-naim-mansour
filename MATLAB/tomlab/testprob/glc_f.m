% glc_f.m
%
% function f = glc_f(x, Prob)
%
% Test functions for constrained global optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Sept 29, 2008.

function f = glc_f(x, Prob)

x=x(:);
P=Prob.P;

if P == 1	% Gomez 2
    f = 0.1*x(1)^2+0.1*x(2)^2;
elseif P == 2	% Gomez 3
    f = (4-2.1*x(1)^2+x(1)^4/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;
elseif P == 3	% HS 59
    u = Prob.user.u;
    f = -u(1)+u(2)*x(1)+u(3)*x(1)^3-u(4)*x(1)^4+u(5)*x(2)-u(6)*x(1)*x(2) ...
        +u(7)*x(2)*x(1)^2+u(8)*x(1)^4*x(2)-u(9)*x(2)^2+u(10)*x(2)^3 ...
        -u(11)*x(2)^4+u(12)/(x(2)+1)+u(13)*x(1)^2*x(2)^2 ...
        +u(14)*x(1)^3*x(2)^2-u(15)*x(1)^3*x(2)^3+u(18)*exp(0.0005*x(1)*x(2)) ...
        -u(19)*x(1)^3*x(2)-u(16)*x(1)*x(2)^2+u(17)*x(1)*x(2)^3-0.12694*x(1)^2;
elseif P == 4	% HS 65
    f = (x(1)-x(2))^2+(x(1)+x(2)-10)^2/9+(x(3)-5)^2;
elseif P == 5	% HS 104
    f = 0.4*x(1)^0.67*x(7)^(-0.67)+0.4*x(2)^0.67*x(8)^(-0.67)+10-x(1)-x(2);
elseif P == 6	% HS 105
    y = Prob.user.y;
    ai = x(1)/x(6)*exp(-(y-x(3)).^2/(2*x(6)^2));
    bi = x(2)/x(7)*exp(-(y-x(4)).^2/(2*x(7)^2));
    ci = (1-x(2)-x(1))/x(8)*exp(-(y-x(5)).^2/(2*x(8)^2));
    f = -sum(log((ai+bi+ci)/sqrt(2*pi)));
elseif P == 7	% S 234
    f = (x(2)-x(1))^4-(1-x(1));
elseif P==8 | P==9 | P==10 | P==23	% S 236, S 237, S 239, S 237 altered
    B = Prob.user.B;
    f = B(1)+B(2)*x(1)+B(3)*x(1)^2+B(4)*x(1)^3+B(5)*x(1)^4+B(6)*x(2)+ ...
        B(7)*x(1)*x(2)+B(8)*x(1)^2*x(2)+B(9)*x(1)^3*x(2)+B(10)*x(1)^4*x(2)+ ...
        B(11)*x(2)^2+B(12)*x(2)^3+B(13)*x(2)^4+B(14)*(1/(x(2)+1))+ ...
        B(15)*x(1)^2*x(2)^2+B(16)*x(1)^3*x(2)^2+B(17)*x(1)^3*x(2)^3+ ...
        B(18)*x(1)*x(2)^2+B(19)*x(1)*x(2)^3+B(20)*(exp(5E-4*x(1)*x(2)));
    f=-f;
elseif P == 11	% S 330
    f = 0.044*x(1)^3/x(2)^2+1/x(1)+0.0592*x(1)/x(2)^3;
elseif P == 12	% S 332
    t = Prob.user.tmp2;
    f = pi/3.6*sum((log(t)+x(2)*sin(t)+x(1)*cos(t)).^2+(log(t)+x(2)*cos(t)-x(1)*sin(t)).^2);
elseif P == 13	% S 343
    f = -0.0201*x(1)^4*x(2)*x(3)^2*1E-7;
elseif P == 14	% FP 3.2 TP 1
    f = x(1)+x(2)+x(3);
elseif P == 15	% FP 3.3 TP 2
    f = 37.293239*x(1)+0.8356891*x(1)*x(5)+5.3578547*x(3)^2-40792.141;
elseif P == 16	% FP 3.5 TP 4
    f = -2*x(1)+x(2)-x(3);
elseif P == 17	% FP 4.10 TP 9
    f = -x(1)-x(2);
elseif P == 18 % Zimmerman
    f = 9-x(1)-x(2);
elseif P == 19 | P == 20 | P == 21  % Bump
    tmp = (1:Prob.N)';
    f = -abs((sum((cos(x)).^4)-2*prod((cos(x)).^2))/(sqrt(tmp'*(x.^2))));
elseif P == 22 % HGO 468:1 + constraint
    f = 4*x(1)*x(2)*sin(4*pi*x(2));
    f=-f;
end

% MODIFICATION LOG
%
% 990205  mbk  First Version
% 990413  mbk  Moved problem 1 to glb_*, added Gomez 2
% 990416  mbk  Small changes in comments.
% 990427  mbk  Problem 'Zimmerman' added.
% 020323  hkh  Adding Floudas-Pardalos chapter 12 problems
% 031126  hkh  Correcting problem 20
% 080603  med  Switched to glcAssign, cleaned
% 080925  nhq  Removed all IP problems to newly created glcIP_c.
