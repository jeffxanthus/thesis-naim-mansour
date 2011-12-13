% jogo09_f.m
%
% function f = jogo09_f(x, Prob)
%
% Test functions for constrained global optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 9, 2009.   Last modified June 10, 2009.

function f = jogo09_f(x, Prob)

x=x(:);
n = length(x);
P=Prob.P;

if P == 1	% Low-dimensional 1 - Mladenovic
    if x(1) < 2
       f = x(2)+10e-5*(x(2)-x(1))^2;
    elseif x(1) < 4
       f = 1/(27*sqrt(3))*((x(1)-3)^2-9)*x(2)^3;
    else
       f = 1/3*(x(1)-2)^3+x(2)-11/3;
    end
elseif P == 2	% Low-dimensional 2 - Sun
    f = -x(1)^2+1/6*x(2)^3-x(1)-43*x(2);
elseif P == 3	% Low-dimensional 3 - Churilov
    f = 1.45*x(1)^2+0.75*x(2)^2+7*x(3)^2-2*x(1)*x(3)-2*x(2)*x(3)+2*sqrt(x(1)+x(2)+x(3));
elseif P == 4	% Box constrained with boundary solutions 1 - Cox
    e = 0.2+0.2*rand(size(x));
    f = sum(2.2*(x+e).^2-(x+e).^4);
elseif P == 5	% Box constrained with boundary solutions 2 - Maranas - Circle Packing
    p = n/2;
    z = [];
    for j = 2:p
       for i = 1:(j-1)
          z = [z;(x(i)-x(j))^2+(x(p+i)-x(p+j))^2];
       end
    end
    f = -max(z);
elseif P == 6	% Box constrained with boundary solutions 3 - Locatelli 1
    f = -3*sum(x.^2)+2*sum(x(1:n-1).*x(2:n));
elseif P == 7	% Box constrained with boundary solutions 3 - Locatelli 2
    f = -sum(1./(1:n)'.*x)*log(1+sum(1./(1:n)'.*x));
elseif P==8 % Box constrained with boundary solutions 3 - Locatelli 3
    f = -sum(x.^2)*log(1+sum(x.^2));
elseif P == 9 % Fully constrained with boundary solutions 1
    f = 5*sum(x(1:4))-5*sum(x(1:4).^2)-sum(x(5:13));
elseif P == 10 | P == 11 | P == 12 % Fully constrained with boundary solutions 2, 3 and 4
    f = x'*Prob.user.Q*x;
elseif P == 13	% Fully constrained with boundary solutions 5
    f = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2 ...
       +4*(x(4)-5)^2+(x(5)-3)^2+2*(x(6)-1)^2+5*x(7)^2 ... % +4*(x(4)-5)^2+(x(3)-3)^2+2*(x(6)-1)^2+5*x(7)^2 in paper, seems to be typo
       +7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+45;
elseif P == 14	% Fully constrained with boundary solutions 6
    f = -prod(sqrt(n)*x);
elseif P == 15	% Fully constrained with boundary solutions 7
    f = -prod(sqrt(n*x));
elseif P == 16	% Fully constrained with boundary solutions 8 - Random constraints
    f = -(x(1)+sum(((2:n)'-1)./(2:n)'.*x(2:n)))^(3/2);
elseif P == 17	% Fully constrained with boundary solutions 9 - Random constraints
    f = -sum(x./(1:n)')*log(1+sum(x./(1:n)'));
elseif P == 18 % Fully constrained with boundary solutions 10 - Random constraints
    f = -3*sum(x.^2)+2*sum(x(1:n-1).*x(2:n));
elseif P == 19 % Fully constrained with boundary solutions 11 - Random constraints
    f = -sum(x.^2)*log(1+sum(x.^2));
elseif P == 20 % Fully constrained with boundary solutions 12 - Random constraints
    f = -sqrt(sum(x.^2))-sqrt(sum((x-1).^2));
elseif P == 21 % Fully constrained with boundary solutions 13 - Random constraints
    f = -sum(log(1+x)-exp(-x/n));
elseif P == 22 % Fully constrained with boundary solutions 14 - Random constraints
    f = -exp(sum(x./(1:n)'))-sqrt(1+(sum(x./(1:n)'))^2);
elseif P == 23 % Fully constrained with boundary solutions 15 - Random constraints
    f = -sqrt(sum(x.^2))-(sum(sqrt(x)))^(3/2);
elseif P == 24 % Fully constrained with boundary solutions 16 - Random constraints and objective
    f = sum(Prob.user.beta.*x.*x(Prob.user.j1).^2.*x(Prob.user.j2).^3);
elseif P == 25 % Fully constrained with boundary solutions 17 - Random constraints and objective
    f = -prod(x.^Prob.user.alfa);
end

% MODIFICATION LOG
%
% 090609  bjo  First version created