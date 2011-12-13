% function f = uc_f(x, Prob)
%
% Test functions for Unconstrained optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Aug 4, 2008.

function f = uc_f(x, Prob)

x=x(:);
P=Prob.P;
uP=Prob.uP;

if P == 1 % Rosenbrock banana
    if isfield(Prob.PartSep,'index')
        i=Prob.PartSep.index;
    else
        i=0;
        Prob.PartSep.pSepFunc=0;
    end
    if i==1 & Prob.PartSep.pSepFunc==2
        % Compute only first term
        f = 100*(x(2)-x(1)^2)^2;
    elseif i==2 & Prob.PartSep.pSepFunc==2
        % Compute only second term
        f = (1-x(1))^2;
    else
        f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    end
elseif P == 2 % Exponential
    f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
elseif P == 3
    v = x(1)*x(2)*x(3);
    f = 50*x(1)*x(2) + 100*x(1)*x(3) + 200*x(2)*x(3) + 25*x(1);
    if v~=0
        f=f+200/v;
    else
        f=f+1E200;
    end
elseif P == 4 % Spiral
    r = norm(x);
    if r==0, r=1E-10; end
    if real(x(2)) >= 0
        phi = acos(x(1)/r);
    else
        phi = 2*pi - acos(x(1)/r);
    end
    f = uP(1)*sin(phi + r) + uP(2)*r;
elseif P == 5 % Fletcher exercise 2.3
    % X=(0,0) saddle; x= 0,-1) saddle; x=(1,0) local min; x =(-1,-1) local max
    f = 2 * x(1)^3 - 3 * x(1)^2 - 6 * x(1) * x(2) * (x(1)-x(2)-1);
elseif P == 6	% Fletcher exercise 3.9
    % Min (2,1). Start (0,0). Conv 2 steps BFGS with exact line search.
    f = x(1)^2+4*x(2)^2-4*x(1)-8*x(2);
elseif P == 7	% AG-1, Anders Goeran 1.
    % 2 min at f(0.354667,-0.354667)=f(-0.354667,0.354667)=1.42052
    f = exp(x(1)^2+x(2)^2)+exp(2*x(1)*x(2)-1)+2*x(1)*x(2)+0.1;
elseif P == 8	% AG-2, Anders Goeran 2.
    % min at f(0,0)=-1. Deep hole, flat function away from origo.
    f = -exp(-x(1)^2-x(2)^2);
elseif P == 9 % RE-1, Roger Enblom, May-96
    % min at f(0,0)=0, max at(0,+-1)=3. Marsden - Tromba fig. 2.1.17
    eVal=exp(1-x(1)^2-x(2)^2);
    f = uP(1)*(x(1)^2+3*x(2)^2)*eVal;
elseif P == 10 % RE-2, Roger Enblom, May-96
    % netlib/trig, n=2, several max, min and saddle points.
    f = 4-sin(x(2))-cos(x(1))-3*cos(x(2));
elseif P == 11	% Structured exponential, p=4
    % total value in position 1, element function values in  positions 2 to p+1
    % Used in testing of STrustR.m
    f = exp(x(1))*[4*x(1)^2; 2*x(2)^2; 4*x(1)*x(2); 2*x(2) + 1];
    if uP(1) == 0
        f=sum(f);
    else
        f = f(uP(1));
    end
elseif P == 12	% Fletcher Q.2.4
    f=(x(2)-x(1)^2)^2+x(1)^5;
elseif P == 13	% Fletcher Q.2.5
    f=2*x(1)^2+x(2)^2-2*x(1)*x(2)+2*x(1)^3+x(1)^4;
elseif P == 14	% Fletcher Q.2.2
    % min in f(0.69,-1.37)=-0.58
    f = x(1)^4+x(1)*x(2)+(1+x(2))^2;
elseif P == 15	% Fletcher Q.2.7 Quadratic
    % 0.5*x'*A*x + b' * x
    f = x(1)^2+4*x(2)^2-4*x(1)-8*x(2);
elseif P == 16	% Fletcher Q.2.6
    f = x(1)^2*x(2)^2-4*x(1)^2*x(2)+4*x(1)^2+2*x(1)*x(2)^2+x(2)^2 ...
        -8*x(1)*x(2)+8*x(1)-4*x(2);
elseif P == 17	% Fletcher Q.3.3
    f = 0.5*(x(1)^2+x(2)^2)*exp(x(1)^2-x(2)^2);
elseif P == 18  % Nash-Sofer page 307 #3
    f=5*x(1)^4 + 6*x(2)^4-6*x(1)^2+2*x(1)*x(2)+5*x(2)^2+15*x(1)-7*x(2)+13;
elseif P == 19  % 3rd degree
    f=5*x(1)^2 +x(2)^2-2*x(1)^3-x(1)*x(2)^2;
end

% MODIFICATION LOG
%
% 981021  hkh  Added handling of Partial Separable functions to P==1
% 981102  hkh  Check if pSepFunc == 2
% 080603  med  Switched to conAssign, cleaned
