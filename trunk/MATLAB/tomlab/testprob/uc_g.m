% function g = uc_g(x,Prob)
%
% Gradient for test functions for Unconstrained optimization.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Aug 4, 2008.

function g = uc_g(x,Prob)

x=x(:);
P=Prob.P;
uP=Prob.uP;

if P == 1
    if isfield(Prob.PartSep,'index')
        i=Prob.PartSep.index;
    else
        i=0;
        Prob.PartSep.pSepFunc=0;
    end
    if i==1 & Prob.PartSep.pSepFunc==2
        % Compute only first term
        g = [-400*x(1)*(x(2)-x(1)^2); 200*(x(2)-x(1)^2)	];
    elseif i==2 & Prob.PartSep.pSepFunc==2
        % Compute only second term
        g = [-2*(1-x(1)); 0];
    else
        g = [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 200*(x(2)-x(1)^2)	];
    end
elseif P == 2
    e=exp(x(1));
    g = [e*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 8*x(1) + 6*x(2) + 1);
        e*(4*x(1) + 4*x(2) + 2)	];
elseif P == 3
    v = x(1)*x(2)*x(3);
    if real(v)~=0
        g = [50*x(2) + 100*x(3) + 25 - 200/(x(1) * v)
            50*x(1) + 200*x(3)      - 200/(x(2) * v)
            100*x(1) + 200*x(2)      - 200/(x(3) * v)];
    else
        g = [50*x(2) + 100*x(3) + 25 - 1000
            50*x(1) + 200*x(3)      - 1000
            100*x(1) + 200*x(2)      - 1000 ];
    end
elseif P == 4
    r = norm(real(x));
    if r==0, r=1E-10; end
    if real(x(2)) >= 0
        phi = acos(x(1)/r);
    else
        phi = 2*pi - acos(x(1)/r);
    end
    c=x(1)/r;
    s=x(2)/r;
    if real(s)==0
        g = uP(1)*cos(phi + r)*[c ; 0] + uP(2)*[c; s];
    else
        g=uP(1)*cos(phi + r)*[c - abs(s)/r; s*(1 + c/(abs(s)*r))] + uP(2)*[c; s];
    end
elseif P == 5
    g = [x(1)^2 - x(1) - 2 * x(1) * x(2) + x(2)*(x(2)+1) ;...
        -x(1)^2 + x(1) + 2 * x(1) * x(2) ];
    g=6*g;
elseif P == 6
    g = [2*x(1)-4;8*x(2)-8];
elseif P == 7	% AG-1, Anders Goeran 1.
    % 2 min at f(0.354667,-0.354667)=f(-0.354667,0.354667)=1.42052
    e1 = exp(x(1)^2+x(2)^2);
    e2 = exp(2*x(1)*x(2)-1);
    g = [2*x(1)*e1+2*x(2)*e2+2*x(2); 2*x(2)*e1+2*x(1)*e2+2*x(1)];
elseif P == 8	% AG-2, Anders Goeran 2.
    % min at f(0,0)=-1. Deep hole, flat function away from origo.
    e1= exp(-x(1)^2-x(2)^2);
    g = [2*x(1)*e1; 2*x(2)*e1];
elseif P == 9 % RE-1
    %Marsden - Tromba fig. 2.1.17min at f(0,0)=0, max at(0,+-1)=3.
    eVal=exp(1-x(1)^2-x(2)^2);
    g = uP(1)*eVal*[2*x(1)*(1-x(1)^2-3*x(2)^2); 2*x(2)*(3-x(1)^2-3*x(2)^2)];
elseif P == 10% RE-2
    %netlib/trig, n=2, several max, min and saddle points.
    g = [sin(x(1)); 3*sin(x(2))-cos(x(2))];
elseif P == 11% Structured exponential, p=4
    eX=exp(x(1));
    if uP(1) ==0
        g = eX*[4*x(1)^2 + 8*x(1) + 2*x(2)^2 + 4*x(1)*x(2) + 6*x(2) + 1; ...
            4*x(2)+4*x(1)+2];
    elseif uP(1)==1
        g = [eX*(4*x(1)^2 + 8*x(1));0 ];
    elseif uP(1)==2
        g = eX*[2*x(2)^2;4*x(2) ];
    elseif uP(1)==3
        g = eX*[4*x(1)*x(2) + 4*x(2); 4*x(1) ];
    elseif uP(1)==4
        g = eX*[2*x(2) + 1 ;2];
    end
elseif P == 12	% Fletcher Q.2.4
    g = [ 4*x(1)^3+5*x(1)^4-4*x(1)*x(2) ; 2*x(2)-2*x(1)^2 ];
elseif P == 13	% Fletcher Q.2.5
    g = [ 4*x(1)-2*x(2)+6*x(1)^2+4*x(1)^3 ; 2*x(2)-2*x(1) ];
elseif P == 14	% Fletcher Q.2.2
    % min i f(0.69,-1.37)=-0.58
    g = [4*x(1)^3+x(2);x(1)+2+2*x(2)];
elseif P == 15	% Fletcher Q.2.7 Quadratic
    g = [2*x(1)-4;8*x(2)-8];
elseif P == 16	% Fletcher Q.2.6
    g = [2*x(1)*x(2)^2-8*x(1)*x(2)+8*x(1)+2*x(2)^2-8*x(2)+8;...
        2*x(1)^2*x(2)-4*x(1)^2+4*x(1)*x(2)+2*x(2)-8*x(1)-4];
elseif P == 17	% Fletcher Q.3.3
    e = exp(x(1)^2-x(2)^2);
    g = e*[x(1)*(1+x(1)^2+x(2)^2); x(2)*(1-x(1)^2-x(2)^2)];
elseif P == 18  % Nash-Sofer page 307 #3
    g=[20*x(1)^3-12*x(1)+2*x(2)+15;24*x(2)^3+2*x(1)+10*x(2)-7];
elseif P == 19  % 
    g=[10*x(1)-6*x(1)^2-x(2)^2;2*x(2)-2*x(1)*x(2)];
end

% MODIFICATION LOG
%
% 981021  hkh  Added handling of Partial Separable functions to P==1
% 981102  hkh  Check if pSepFunc == 2
% 080603  med  Switched to conAssign, cleaned
