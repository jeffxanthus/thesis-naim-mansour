% Defines matrix of constraint gradients for mixed-integer
% nonlinear programming (MINLP) problems.
%
% function dc = minlp_dc(x,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function dc = minlp_dc(x,Prob)

x=x(:);
P = Prob.P;

if P==1
    dc = [ ...
        2*x(1)     0.0       1.0  0.0  0.0 ; ...
        0.0   1.5*sqrt(x(2)) 0.0  1.5  0.0];
elseif P==10
    A = [1 0     1  0 0 ; ...
        0 1.333 0  1 0 ; ...
        0 0    -1 -1 1 ];
    dc = [ ...
        2*x(1)     0.0       1.0  0.0  0.0 ; ...
        0.0   1.5*sqrt(x(2)) 0.0  1.5  0.0;A];
elseif P==2
    dc = [-exp(x(1)-0.2),-1,0];
elseif P==3
    y = x(1:4);
    x = x(5:7);
    dc = [0    0    2*y(3)   0   2*x(1) 2*x(2) 2*x(3) ; ...
        0 2*y(2)     0     0     0    2*x(2)   0    ; ...
        0    0    2*y(3)   0     0      0    2*x(3) ; ...
        0 2*y(2)     0     0     0      0    2*x(3) ];
    % If ConsPattern is given, dc must either be sparse, or define only
    % the non-zero values in a vector
elseif P==9
    y = x(1:4);
    x = x(5:7);
    dc = [0    0    2*y(3)   0   2*x(1) 2*x(2) 2*x(3) ; ...
        0 2*y(2)     0     0     0    2*x(2)   0    ; ...
        0    0    2*y(3)   0     0      0    2*x(3) ; ...
        0 2*y(2)     0     0     0      0    2*x(3) ];
    % Linear constraints treated as nonlinear
    dc = [dc; ...
        1 1 1 0 1 1 1 ; ...
        1 0 0 0 1 0 0 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 0 0 1 1 0 0 ];
elseif P==4
    dc = [-7+1.2*x(1)^0.2*x(2)^1.7 , -9+1.7*x(1)^1.2*x(2)^0.7 ];
elseif P>=5 & P<=8
    p = Prob.user.P;
    N = Prob.user.N;
    m = x(1:p);
    r = reshape(x(2*p+1:end),N,p);
    dc = [r' zeros(N,p)];
    for k=1:p
        dc(k,(k-1)*p+(2*p+1:3*p))=m';
    end
elseif P==11
    z11 = 0.9*x(2)*0.5*exp(-0.5*x(1));
    z12 = 0.9*(1-exp(-0.5*x(1)));
    z21 = 0.8*x(2)*0.4*exp(-0.4*x(1 ));
    z22 = 0.8*(1-exp(-0.4*x(1)));
    M  = Prob.user.M;
    % Use big-M-method and equalities rewritten as two inequalities
    dc = [ z11,  z12, M, 0 ; ...
        -z11, -z12, M, 0 ; ...
        z21,  z22, 0, M ; ...
        -z21, -z22, 0, M];
elseif P==12
    z11 = 2*x(1);
    z12 = 2*x(2);
    z21 = 2*(x(1)-4);
    z22 = 2*(x(2)-1);
    z31 = 2*(x(1)-2);
    z32 = 2*(x(2)-4);
    M  = Prob.user.M;
    % Use big-M-method
    dc = [ z11,  z12, M, 0, 0 ; ...
        z21,  z22, 0, M, 0 ; ...
        z31,  z32, 0, 0, M];
elseif P==13
    dc = [log(0.1) log(0.2) log(0.15) 0 0 0 0 0 1/(1-x(9)) 0 0;...
        0 0 0 log(0.05) log(0.2) log(0.15) 0 0 0 1/(1-x(10)) 0;...
        0 0 0 0 0 0 log(0.02) log(0.06) 0 0 1/(1-x(11))];
elseif P==14
    dc = [0 0 0 -x(4)^(-0.5)*x(5)^2+8  4*x(5)-x(5)^(-0.5)-4*x(4)^0.5*x(5)+11];
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 080603 med  Switched to minlpAssign, cleaned