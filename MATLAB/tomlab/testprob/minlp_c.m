% Defines constraints for mixed-integer nonlinear programming (MINLP) problems.
%
% function c = minlp_c(x,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function c = minlp_c(x,Prob)

x=x(:);
P = Prob.P;

if P==1
    c = [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4)];
elseif P==10
    A = [1 0     1  0 0 ; ...
        0 1.333 0  1 0 ; ...
        0 0    -1 -1 1 ];
    c = [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4);A*x];
elseif P==2
    c = -exp(x(1)-0.2)-x(2);
elseif P==3
    y = x(1:4);
    x = x(5:7);
    c = [ y(3)^2+sum(x.^2) ; ...
        y(2)^2+x(2)^2 ; ...
        y(3)^2+x(3)^2 ; ...
        y(2)^2+x(3)^2 ];
elseif P==9
    % Linear c/s treated as general nonlinear
    A = [ 1 1 1 0 1 1 1 ; ...
        1 0 0 0 1 0 0 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 0 0 1 1 0 0 ];
    c = A*x;
    y = x(1:4);
    x = x(5:7);
    c = [ y(3)^2+sum(x.^2) ; ...
        y(2)^2+x(2)^2 ; ...
        y(3)^2+x(3)^2 ; ...
        y(2)^2+x(3)^2;c ];
elseif P==4
    c = (x(1)^1.2)*(x(2)^1.7) - 7*x(1) - 9*x(2);
elseif P>=5 & P<=8
    p = Prob.user.P;
    N = Prob.user.N;
    m = x(1:p);
    r = reshape( x(2*p+1:end), N, p);
    c = r'*m;
elseif P==11
    z1 = 0.9*x(2)*(1-exp(-0.5*x(1))) - 10;
    z2 = 0.8*x(2)*(1-exp(-0.4*x(1))) - 10;
    M  = Prob.user.M;
    % Use big-M-method and equalities rewritten as two inequalities
    c = [ z1-M*(1-x(3)) ; ...
        -z1-M*(1-x(3)) ; ...
        z2-M*(1-x(4)) ; ...
        -z2-M*(1-x(4)) ...
        ];
    c(isnan(c)) = 1E3;
elseif P==12
    z1 = x(1)^2 + x(2)^2 -1;
    z2 = (x(1)-4)^2 + (x(2)-1)^2 -1;
    z3 = (x(1)-2)^2 + (x(2)-4)^2 -1;
    M  = Prob.user.M;
    % Use big-M-method
    c = [ z1-M*(1-x(3)) ; ...
        z2-M*(1-x(4)) ; ...
        z3-M*(1-x(5))];
elseif P==13
    c(3,1) = -log(1-x(11)) + log(0.02)*x(7) + log(0.06)*x(8);
    c(1,1) = -log(1-x(9))  + log(0.1)*x(1)  + log(0.2)*x(2) + log(0.15)*x(3);
    c(2,1) = -log(1-x(10)) + log(0.05)*x(4) + log(0.2)*x(5) + log(0.15)*x(6);
elseif P==14
    c = 2*x(5)^2-2*x(5)^0.5-2*x(4)^0.5*x(5)^2+11*x(5)+8*x(4);
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 070307 hkh  Problem 11 may generate NaN values for certain x
% 080603 med  Switched to minlpAssign, cleaned