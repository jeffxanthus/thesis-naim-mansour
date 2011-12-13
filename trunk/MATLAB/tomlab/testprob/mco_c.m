% function cx=mco_c(x, Prob)
%
% The mc problems are Uncontrained & Constrained Nonlinear Multi Criterium Problems
%
% mco_c evaluates the nonlinear constraints for test problem P (Prob.P)
% at the point x.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function cx=mco_c(x, Prob)

P=Prob.P;
cx=[];

if P==1
    % 'MCO-TP 1'
elseif P==2
    % 'MCO-TP 2'
elseif P==3
    % 'MCO-TP 3'
elseif P==4
    % 'MCO-TP 4'
    p3=(x(1)-1)^2+(x(2)+0.5)^2;
    cx = [-(1-exp(-p3)); 1-exp(-p3)];
elseif P==5
    % 'MCO-TP 5'
elseif P==6
    % 'MCO-TP 6'
elseif P==7
    % 'MCO-TP 7'
elseif P==8
    % 'MCO-TP 8'
    cx = x(1)^2+x(2)^2;
elseif P==9
    % 'MCO-TP 9'
    cx = [-(x(3)-3)^2-x(4);
        (x(5)-3)^2+x(6)];
end

% MODIFICATION LOG:
%
% 040511  med  Created
% 080603  med  Switched to clsAssign, cleaned