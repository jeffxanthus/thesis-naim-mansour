% glc_c.m
%
% function [cx]=glc_c(x, Prob)
%
% glc_c evaluates the constraints for test problem Prob.P
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Sept 29, 2008.

function cx = glc_c(x, Prob)

P=Prob.P;

switch P
    case 1 % Gomez 2
        cx = -sin(4*pi*x(1))+2*(sin(2*pi*x(2)))^2;
    case 2 % Gomez 3
        cx = -sin(4*pi*x(1))+2*(sin(2*pi*x(2)))^2;
    case 3 % HS 59
        cx = [x(1)*x(2)-700;x(2)-x(1)^2/125;(x(2)-50)^2-5*(x(1)-55)];
    case 4 % HS 65
        cx = 48-(x(1)^2+x(2)^2+x(3)^2);
    case 5 % HS 104
        cx = [1-0.0588*x(5)*x(7)-0.1*x(1)
            1-0.0588*x(6)*x(8)-0.1*(x(1)+x(2))
            1-4*x(3)/x(5)-2*x(3)^(-0.71)/x(5)-0.0588*x(3)^(-1.3)*x(7)
            1-4*x(4)/x(6)-2*x(4)^(-0.71)/x(6)-0.0588*x(4)^(-1.3)*x(8)
            0.4*x(1)^0.67*x(7)^(-0.67)+0.4*x(2)^0.67*x(8)^(-0.67)+10-x(1)-x(2)
            -(0.4*x(1)^0.67*x(7)^(-0.67)+0.4*x(2)^0.67*x(8)^(-0.67)+10-x(1)-x(2))];
    case 6 % HS 105
        cx = [];
    case 7 % TP 234
        cx = -x(1)^2-x(2)^2+1;
    case 8 % TP 236
        cx = [x(1)*x(2)-700;x(2)-5*(x(1)/25)^2];
    case 9 % TP 237
        cx = [x(1)*x(2)-700;x(2)-5*(x(1)/25)^2;(x(2)-50)^2-5*(x(1)-55)];
    case 10 % TP 239
        cx = x(1)*x(2)-700;
    case 11 % TP 330
        cx = 1-8.62*x(2)^3/x(1);
    case 12 % TP 332
        t = Prob.user.tmp2;
        pmax = max(180/pi*atan(abs((1./t-x(1))./(log(t)+x(2)))));
        cx = pmax;
    case 13 % TP 343
        cx = [675-x(1)^2*x(2);0.419-1E-7*x(1)^2*x(3)^2];
    case 14 % FP 3.2  - TP 1
        cx = [100*x(1)-x(1)*x(6)+833.33252*x(4)-83333.333
            x(2)*x(4)-x(2)*x(7)-1250*x(4)+1250*x(5)
            x(3)*x(5)-x(3)*x(8)-2500*x(5)];
    case 15 % FP 3.3  - TP 2
        cx = [-0.0022053*x(3)*x(5)+0.0056858*x(2)*x(5)+0.0006262*x(1)*x(4)
            0.0071317*x(2)*x(5)+0.0021813*x(3)^2+0.0029955*x(1)*x(2)
            0.0047026*x(3)*x(5)+0.0019085*x(3)*x(4)+0.0012547*x(1)*x(3)];
    case 16 % FP 3.5  - TP 4
        B=[0 0 1;0 -1 0;-2 1 -1];
        r=[1.5 -0.5 -5]';
        cx = x'*B'*B*x-2*r'*B*x;
    case 17 % FP 4.10 - TP 9
        cx = [2*x(1)^4-8*x(1)^3+8*x(1)^2-x(2);4*x(1)^4-32*x(1)^3+88*x(1)^2-96*x(1)-x(2)];
    case 18 % Zimmermann
        cx = [(x(1)-3)^2+(x(2)-2)^2;x(1)*x(2)];
    case {19,20,21} % Bump
        cx = prod(x);
    case 22 % HGO 468:1 + constraint
        cx = x(1)-0.7-15*(x(2)-0.6)^2;
        %cx = x(1)-0.7-0.5*(x(2)-0.5)^2;
        %cx = x(1)^2-x(2);
    case 23 % TP 237 altered
        cx = [(x(2)-75)^2-5*(x(1)-55); x(2)+0.2*(x(1)-35).^2-70];
end

% MODIFICATION LOG
%
% 990413  mbk  Moved problem 1 to glb_*, added Gomez 2
% 990416  mbk  Small changes in comments.
% 990427  mbk  Problem 'Zimmerman' added.
% 020323  hkh  Adding Floudas-Pardalos chapter 12 problems
% 080603  med  Switched to glcAssign, cleaned
% 080925  nhq  Removed all IP problems to newly created glcIP_c.
