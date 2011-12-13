% function f = lgo1_f(x, Prob)
%
% Test functions for global optimization. One dimension.

% Reference:
% Pintér, J.D., Bagirov, A., and Zhang, J. (2003) An Illustrated Collection of
% Global Optimization Test Problems. Research Report, Pintér Consulting Services,
% Inc. Halifax, NS, Canada; and CIAO-ITMS, University of Ballarat, Ballarat, Vic.,
% Australia.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function f = lgo1_f(x, Prob)

x=x(:);
P=Prob.P;

if P == 1 % P1
    f = x.^4 - 3*x.^3 - 8*x.^2 - 20*x + 3;
elseif P == 2 % P2
    f = 0.2*x.^5 - 1.6995*x.^4 + 0.998266*x.^3 - 0.0218343*x.^2 - 0.000089248*x;
elseif P == 3 % P3
    f = 0.1666667*x.^6 - 1.05*x.^4 + 2*x.^2 - x;
elseif P == 4 % P4
    f = 24*x.^4 - 142*x.^3 + 303*x.^2 - 276*x + 3;
elseif P == 5 % P5
    f = ((x - 1)*(x - 3)*(x - 5)*(x - 7)).^2;
elseif P == 6 % T1
    f = sin(x)*(cos(x)-sin(x)).^2;
elseif P == 7 % T2
    f = sin(x.^2+x)+cos(3*x);
elseif P == 8 % T3
    f=0;
    for i=1:5
        f = f + i*sin((i+1)*x+i);
    end
elseif P == 9 % T4
    f = exp(-x)-sin(x).^3;
elseif P == 10 % T5
    f=0.1*x+sqrt(abs(x))*sin(x).^2;
elseif P == 11 % T6
    f=x.^2-cos(18*x);
elseif P == 12 % T7
    f=-20*exp(-0.02*abs(x))-exp(cos(2*pi*x))+20+exp(1);
elseif P == 13 % T8
    f=x.^2-0.1*cos(5*pi*x);
elseif P == 14 % T9
    f=-sin(x)*(sin(x.^2/pi)).^3;
elseif P == 15 % T10
    f=1+x.^2/4000-cos(x);
elseif P == 16 % T11
    f=0.1*(sin(3*pi*x).^2+(x-1).^2*(1+sin(2*pi*x).^2));
elseif P == 17 % T12
    f=-(0.806*cos((x-9.681).^2/pi)*exp(-pi*(x-9.681).^2) + 0.517*cos((x-9.4).^2/pi)*exp(-pi*(x-9.4).^2) + 0.1*cos((x-8.025).^2/pi)*exp(-pi*(x-8.025).^2) + 0.908*cos((x-2.196).^2/pi)*exp(-pi*(x-2.196).^2) + 0.965*cos((x-8.074).^2/pi)*exp(-pi*(x-8.074).^2));
elseif P == 18 % T13
    f=0;
    for i=1:3
        f = f + (-1).^i*sin(2*pi*i*x);
    end
elseif P == 19 % T14
    f=exp(-3*x)-sin(x).^3;
elseif P == 20 % T15
    f=0;
    for i=1:5
        f = f + (4-cos((i+1)*x));
    end
elseif P == 21 % T16
    f=(3*x-1.4)*sin(18*x)+1.7;
elseif P == 22 % T17
    f=cos(x)-sin(5*x)+1;
elseif P == 23 % T18
    f=-x-sin(3*x)+1.6;
elseif P == 24 % T19
    f=sin(x)*x-.84*x+log(x)+sin(10*x/3)+1.3;
elseif P == 25 % T20
    f=0;
    for i=0:5
        f = f + i*cos(i+(i+1)*x)+12;
    end
elseif P == 26 % T21
    f=0.03*(x-10).^2-sin(x-30)+sin(2*(x-20));
elseif P == 27 % T22
    f=sin(3*x)+sin(5*x)+sin(7*x)+sin(11*x);
end

% MODIFICATION LOG
%
% 040114  med  Created
% 080603  med  Switched to glcAssign, cleaned