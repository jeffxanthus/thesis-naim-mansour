% glcIP_f.m
%
% function f = glcIP_f(x, Prob)
%
% Test functions for constrained global optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2008-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Sept 25, 2008.   Last modified Sept 29, 2008.

function f = glcIP_f(x, Prob)

x=x(:);
P=Prob.P;

if P == 1	    % FP 12.2 TP 1
    y = x(1:3); % Integer variables
    f = 2*x(4)+3*x(5)+1.5*y(1)+2*y(2)-0.5*y(3);
elseif P == 2	% FP 12.2 TP 2
    y = x(1); % Integer variables
    f = -0.7*y + 5*(x(2)-0.5)^2 + 0.8;
elseif P == 3	% FP 12.2 TP 3
    y = x(4:7); % Integer variables
    f = (y(1)-1)^2 + (y(2)-2)^2 + (y(3)-1)^2 - log(y(4)+1) + ...
        (x(1)-1)^2 + (x(2)-2)^2 + (x(3)-3)^2;
elseif P == 4	% FP 12.2 TP 4
    f = -x(1)*x(2)*x(3);
elseif P == 5   % FP 12.2 TP 5
    f = 7*x(1) + 10*x(2);
elseif P == 6 	% FP 12.2 TP 6
    % Second variable is integer
    f = -5*x(1) + 3*x(2);
elseif P == 7   % FP 3.3 TP 2 IP
    f = 37.293239*x(1)+0.8356891*x(1)*x(5)+5.3578547*x(3)^2-40792.141;
elseif P == 8	% FP 3.4 TP 3
    f = -25*(x(1)-2)^2-(x(2)-2)^2-(x(3)-1)^2-(x(4)-4)^2-(x(5)-1)^2-(x(6)-4)^2;
elseif P == 9	% FP 3.5 TP 4 IP
    f = -2*x(1)+x(2)-x(3);
elseif P == 10  % FP 4.10 TP 9 IP
    f = -x(1)-x(2);
elseif P == 11	% Kocis & Grossmann 1989
    f = 7.5*x(3) + 5.5*x(4) + 7*x(1)*x(3) + 6*x(1)*x(4) + 5*x(2);
elseif P == 12  % Kocis & Grossmann 1998
    f = [2 3 1.5 2 -0.5]*x;
elseif P == 13  % Lee & Grossmann Disjunctive B&B
    f = (x(1) - 3)^2 + (x(2) - 2)^2 + 2*x(3) + x(4) + 3*x(5);
elseif P == 14  % Kesavan et al. 2004 D
    f = -5*x(4)+ 3*x(5);
elseif P == 416	% FP 12.2 TP 4 GAMS version
    f = -x(1)*x(2)*x(3);
elseif P == 417	% FP 12.2 TP 4 GAMS Int 1st
    f = -x(9)*x(10)*x(11);
elseif P == 418	% FP 12.2 TP 4 GAMS Int last
    y = x(3:5); % Integer variables
    f = 2*x(1)+3*x(2)+1.5*y(1)+2*y(2)-0.5*y(3);
end

% MODIFICATION LOG
%
% 080925  nhq  File created. IP-Problems merged from glc_prob and minlp_prob
% 080929  nhq  File re-created. Names preserved from the old files.

