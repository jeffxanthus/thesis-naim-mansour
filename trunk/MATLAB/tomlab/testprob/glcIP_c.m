% glcIP_c.m
%
% function [cx]=glcIP_c(x, Prob)
%
% glcIP_c evaluates the constraints for test problem Prob.P
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2008-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Sept 25, 2008.   Last modified Sept 29, 2008.

function cx = glcIP_c(x, Prob)

P=Prob.P;

switch P
    case 1  % FP 12.2 TP 1
        y  = x(1:3); % Integer variables
        cx = [x(4)^2+y(1); x(5)^1.5 + 1.5* y(2)];
    case 2  % FP 12.2 TP 2
        cx = -exp(x(2)-0.2) - x(3);
    case 3  % FP 12.2 TP 3
        y  = x(4:7); % Integer variables
        cx = [x(1)^2+x(2)^2+x(3)^2 + y(3)^2; ...
            x(2)^2+y(2)^2; x(3)^2+y(3)^2; x(3)^2+y(2)^2];
    case 4  % FP 12.2 TP 4
        y  = x(4:11); % Integer variables
        cx = [x(1) + 0.1^y(1)*0.2^y(2)*0.15^y(3);
              x(2) + 0.05^y(4)*0.2^y(5)*0.15^y(6);
              x(3)+0.02^y(7)*0.06^y(8)];
    case 5  % FP 12.2 TP 5
        cx = x(1)^1.2*x(2)^1.7-7*x(1)-9*x(2);
    case 6 % FP 12.2 TP 6
        y = x(2);
        cx = 2*y^2 - 2*sqrt(y) - 2*sqrt(x(1))*y^2+11*y+8*x(1);
    case 7  % FP 3.3 TP 2 IP
        cx = [-0.0022053*x(3)*x(5)+0.0056858*x(2)*x(5)+0.0006262*x(1)*x(4)
            0.0071317*x(2)*x(5)+0.0021813*x(3)^2+0.0029955*x(1)*x(2)
            0.0047026*x(3)*x(5)+0.0019085*x(3)*x(4)+0.0012547*x(1)*x(3)];
    case 8  % FP 3.4 TP 3
        cx = [(x(3)-3)^2+x(4);(x(5)-3)^2+x(6)];
    case 9  % FP 3.5 TP 4 IP
        B=[0 0 1;0 -1 0;-2 1 -1];
        r=[1.5 -0.5 -5]';
        cx = x'*B'*B*x-2*r'*B*x;
    case 10 % FP 4.10 - TP 9 IP
        cx = [2*x(1)^4-8*x(1)^3+8*x(1)^2-x(2);4*x(1)^4-32*x(1)^3+88*x(1)^2-96*x(1)-x(2)];
    case 11 % Kocis & Grossmann 1989
        z1 = 0.9*x(2)*(1-exp(-0.5*x(1))) - 10;
        z2 = 0.8*x(2)*(1-exp(-0.4*x(1))) - 10;
        M  = Prob.user.M;
        % Use big-M-method and equalities rewritten as two inequalities
        cx = [ z1-M*(1-x(3)) ; ...
              -z1-M*(1-x(3)) ; ...
               z2-M*(1-x(4)) ; ...
              -z2-M*(1-x(4)) ...
             ];
        cx(isnan(cx)) = 1E3;
    case 12 % Kocis & Grossmann 1998
        cx = [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4)];
    case 13 % Lee & Grossmann Disjunctive B&B
        z1 = x(1)^2 + x(2)^2 -1;
        z2 = (x(1)-4)^2 + (x(2)-1)^2 -1;
        z3 = (x(1)-2)^2 + (x(2)-4)^2 -1;
        M  = Prob.user.M;
        % Use big-M-method
        cx = [z1-M*(1-x(3)) ; ...
              z2-M*(1-x(4)) ; ...
              z3-M*(1-x(5))];
    case 14  % Kesavan et al. 2004 D
        cx = 2*x(5)^2-2*x(5)^0.5-2*x(4)^0.5*x(5)^2+11*x(5)+8*x(4);
    case 416 % FP 12.2 TP 4 GAMS version
        y  = x(4:11); % Integer variables
        cx = [-log(1-x(1)) + log(0.1)*y(1)+log(0.2)*y(2)+log(0.15)*y(3); ...
            -log(1-x(2)) + log(0.05)*y(4)+log(0.2)*y(5)+log(0.15)*y(6);  ...
            -log(1-x(3)) + log(0.02)*y(7)+log(0.06)*y(8)];
    case 417 % FP 12.2 TP 4 GAMS Int 1st
        y  = x(1:8); % Integer variables
        cx = [-log(1-x(9))  + log(0.1)*y(1)+log(0.2)*y(2)+log(0.15)*y(3); ...
            -log(1-x(10)) + log(0.05)*y(4)+log(0.2)*y(5)+log(0.15)*y(6);  ...
            -log(1-x(11)) + log(0.02)*y(7)+log(0.06)*y(8)];
    case 418 % FP 12.2 TP 4 GAMS Int last
        y  = x(3:5); % Integer variables
        cx = [x(1)^2+y(1); x(2)^1.5 + 1.5* y(2)];
end

% MODIFICATION LOG
%
% 080925  nhq  File created. IP-Problems merged from glc_prob and minlp_prob
% 080929  nhq  File re-created. Names preserved from the old files.
