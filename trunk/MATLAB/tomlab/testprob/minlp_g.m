% Defines gradient of objective function for mixed-integer nonlinear
% programming (MINLP) problems.
%
% function g = minlp_g(x,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function g = minlp_g(x,Prob)

x=x(:);
P = Prob.P;

if P==1 | P==10
    g = [2 3 1.5 2 -0.5]';
elseif P==2
    g = [10*(x(1)-0.5);0;-0.7];
elseif P==3 | P==9
    g = [ 2*(x(1)-1) ; ...
        2*(x(2)-2) ; ...
        2*(x(3)-1) ; ...
        -1/(x(4)+1) ; ...
        2*(x(5)-1) ; ...
        2*(x(6)-2) ; ...
        2*(x(7)-3)];
elseif P==4
    g = [7;10];
elseif P>=5 & P<=8
    p = Prob.user.P;
    N = Prob.user.N;
    g = [ ...
        ones(p,1)  ; ...
        0.1*(1:p)' ; ...
        zeros(N*p,1) ];
elseif P==11
    g = [7*x(3)+6*x(4); 5 ; 7.5  + 7*x(1); 5.5 + 6*x(1)];
elseif P==12
    g = [2*(x(1) - 3); 2*(x(2) - 2); 2; 1 ; 3];
elseif P==13
    g = [0;0;0;0;0;0;0;0;-x(10)*x(11); -x(9)*x(11); -x(9)*x(10)];
elseif P==14
    g = [0; 0; 0; -5; 3];
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 080603 med  Switched to minlpAssign, cleaned