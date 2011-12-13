% function dc=goals_dc(x, Prob)
%
% The goals problems are
% Multi Criterium Unconstrained & Constrained Nonlinear Problems
%
% goals_dc computes the gradient to the nonlinear constraints c in the
% point x for the test problem P (Prob.P)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function dc=goals_dc(x, Prob)

P=Prob.P;
dc=[];

if P==1
    %'EASY-FIT TP269'
elseif P==2
    dc=[-3*(1-x(1))^2 -1];
elseif P==3
    %'EASY-FIT TP46'
    dc=[2*x(1)*x(4) 0  0  x(1)^2+cos(x(4)-x(5))  -cos(x(4)-x(5));
        0  1  4*(x(3)^3)*x(4)^2  (x(3)^4)*2*x(4) 0];
elseif P==4
    %'EASY-FIT TP48'
elseif P==5
    %'EASY-FIT TP354'
elseif P==6
    %'EASY-FIT TP355'
    dc1=30*x(2)^2*x(4)^2*x(3)+60*x(2)*x(4)*x(1)-12*x(2)*x(4)*x(3)-126*x(2)^2*x(4)^2*x(1)-6*x(1)...
        +60*x(4)-60*x(2)-30*x(1)*x(4)^2+6*x(3)*x(4)^2+300*x(2)^2*x(4);
    dc2=600*x(2)*x(4)*x(1)-6*x(2)*x(4)^2*x(3)^2-12*x(1)*x(3)*x(4)-120*x(2)*x(4)*x(3)-126*x(2)*x(4)^2*x(1)^2 ...
        -60*x(1)+6*x(4)-600*x(2)-30*x(2)*x(4)^2+30*x(1)^2*x(4)+60*x(2)*x(4)^2*x(3)*x(1);
    dc3=-6*x(2)^2*x(4)^2*x(3)-12*x(2)*x(4)*x(1)+30*x(2)^2*x(4)^2*x(1)+1.5*x(3)+6*x(1)*x(4)^2-60*x(2)^2*x(4);
    dc4=-6*x(2)^2*x(4)*x(3)^2+12*x(1)*x(3)*x(4)-126*x(2)^2*x(4)*x(1)^2-12*x(1)*x(2)*x(3)+60*x(1)-6*x(4)+6*x(2)...
        -30*x(2)^2*x(4)-30*x(1)^2*x(4)+60*x(2)^2*x(4)*x(3)*x(1)+30*x(1)^2*x(2)-60*x(2)^2*x(3)+300*x(1)*x(2)^2;
    dc=[dc1 dc2 dc3 dc4];
elseif P==7
    %'EASY-FIT TP372'
    dc=[1  exp(-5*x(3)) -5*x(2)*exp(-5*x(3))  1 0 0 0 0 0;
        1  exp(-3*x(3)) -3*x(2)*exp(-3*x(3))  0 1 0 0 0 0;
        1  exp(-x(3))   -x(2)*exp(-x(3))      0 0 1 0 0 0;
        1  exp(x(3))     x(2)*exp(x(3))       0 0 0 1 0 0;
        1  exp(3*x(3))   3*x(2)*exp(3*x(3))   0 0 0 0 1 0;
        1  exp(5*x(3))   5*x(2)*exp(5*x(3))   0 0 0 0 0 1;
        -1 -exp(-5*x(3))  5*x(2)*exp(-5*x(3))  1 0 0 0 0 0;
        -1 -exp(-3*x(3))  3*x(2)*exp(-3*x(3))  0 1 0 0 0 0;
        -1 -exp(-x(3))      x(2)*exp(-x(3))    0 0 1 0 0 0;
        -1 -exp(x(3))    -x(2)*exp(x(3))       0 0 0 1 0 0;
        -1 -exp(3*x(3))  -3*x(2)*exp(3*x(3))   0 0 0 0 1 0;
        -1 -exp(5*x(3))  -5*x(2)*exp(5*x(3))   0 0 0 0 0 1];
elseif P==8
    %'EASY-FIT TP373'
    dc=[1  exp(-5*x(3)) -5*x(2)*exp(-5*x(3))  1 0 0 0 0 0;
        1  exp(-3*x(3)) -3*x(2)*exp(-3*x(3))  0 1 0 0 0 0;
        1  exp(-x(3))   -x(2)*exp(-x(3))      0 0 1 0 0 0;
        1  exp(x(3))     x(2)*exp(x(3))       0 0 0 1 0 0;
        1  exp(3*x(3))   3*x(2)*exp(3*x(3))   0 0 0 0 1 0;
        1  exp(5*x(3))   5*x(2)*exp(5*x(3))   0 0 0 0 0 1];
end

% MODIFICATION LOG:
%
% 040517  med  Created.
% 080603  med  Switched to clsAssign, cleaned