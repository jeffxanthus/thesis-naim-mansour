% Defines objective function for mixed-integer nonlinear
% programming (MINLP) problems.
%
% function f = minlp_f(x,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function f = minlp_f(x,Prob)

x=x(:);
P = Prob.P;

if P==1 | P==10
    f = [2 3 1.5 2 -0.5]*x;
elseif P==2
    f = -0.7*x(3)+5*(x(1)-0.5)^2+0.8;
elseif P==3 | P==9
    f = (x(1)-1)^2 + (x(2)-2)^2 + (x(3)-1)^2 - log(x(4)+1) + ...
        (x(5)-1)^2 + (x(6)-2)^2 + (x(7)-3)^2;
elseif P==4
    f = 7*x(1)+10*x(2);
elseif P>=5 & P<=8
    % Trim Loss
    p = Prob.user.P;
    m = x(1:p);
    y = x(p+1:2*p);
    f = sum( m + ((1:p)').*y/10 );
elseif P==11
    f = 7.5*x(3) + 5.5*x(4) + 7*x(1)*x(3) + 6*x(1)*x(4) + 5*x(2);
elseif P==12
    f = (x(1) - 3)^2 + (x(2) - 2)^2 + 2*x(3) + x(4) + 3*x(5);
elseif P==13
    f = -x(9)*x(10)*x(11);
elseif P==14
    f = -5*x(4)+ 3*x(5);
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 080603 med  Switched to minlpAssign, cleaned