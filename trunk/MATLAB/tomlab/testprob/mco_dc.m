% function dc=mco_dc(x, Prob)
%
% The mco problems are Unconstrained & Constrained Nonlinear Multi Criterium
% problems
%
% mco_dc computes the gradient to the nonlinear constraints c in the point x for
% the test problem P (Prob.P)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function dc=mco_dc(x, Prob)

P=Prob.P;
dc=[];

if P==1
    % 'MCO-TP 1'
elseif P==2
    % 'MCO-TP 2'
elseif P==3
    % 'MCO-TP 3'
elseif P==4
    % 'MCO-TP 4'
    p3=(x(1)-1)^2+(x(2)+0.5)^2;
    dc = [-2*(x(1)-1)*exp(-p3) -2*(x(2)+0.5)*exp(-p3);
        2*(x(1)-1)*exp(-p3) 2*(x(2)+0.5)*exp(-p3)];
elseif P==5
    % 'MCO-TP 5'
elseif P==6
    % 'MCO-TP 6'
elseif P==7
    % 'MCO-TP 7'
elseif P==8
    % 'MCO-TP 8'
    dc = [2*x(1) 2*x(2)];
elseif P==9
    % 'MCO-TP 9'
    dc = [0 0 -2*(x(3)-3) -1 0 0;
        0 0 0 0 2*(x(5)-3) 1];
end

% MODIFICATION LOG:
%
% 040511  med  Created
% 080603  med  Switched to clsAssign, cleaned