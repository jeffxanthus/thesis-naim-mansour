% Defines Constraint Hessian matrix for mixed-integer
% nonlinear programming (MINLP) problems.
%
% The matrix is computed by multiplying the constraint Hessian for
% each constraint with the corresponding Lagrange multiplier
% supplied as input lam.
%
% function d2c = minlp_d2c(x,lam,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function d2c = minlp_d2c(x,lam,Prob)

x=x(:);
P = Prob.P;

if P==1 | P==10
    d2c      = spalloc(5,5,2);
    d2c(1,1) = lam(1)*2;
    d2c(2,2) = lam(2)*3/(4*sqrt(x(2)));
elseif P==2
    d2c = spalloc(3,3,1);
    d2c(1,1) = lam(1)*(-exp(x(1)-0.2));
elseif P==3 | P==9
    d = [0 lam(2)+lam(4) lam(1)+lam(3) 0 lam(1) lam(1)+lam(2) lam(1)+lam(3)+lam(4)];
    d2c = sparse(diag(d));
elseif P==4
    d2c = zeros(2,2);
    d2c(1,1) = 0.24*x(1)^(-0.8)*x(2)^1.7;
    d2c(1,2) = 2.04*x(1)^0.2*x(2)^0.7;
    d2c(2,1) = d2c(1,2);
    d2c(2,2) = 1.19*x(1)^1.2*x(2)^(-0.3);
    d2c = d2c*lam(1);
elseif P>=5 & P<=8
    p = Prob.user.P;
    N = Prob.user.N;
    % Zero matrix of size pxp
    Zp = spalloc(p,p,0);
    % Unit matrix pxp
    Ip = speye(p,p);
    d2c = spalloc( (2+N)*p, (2+N)*p, 0);
    for k=1:length(lam)
        % Block-row of zeros with a unit matrix in block k
        Ik = [repmat(Zp,1,k-1) Ip repmat(Zp,1,p-k) ];
        d2c = d2c + lam(k)* ...
            [ repmat(Zp,1,2)   Ik; ...
            spalloc(p,(2+N)*p,0) ; ... % entire block-row of zeros
            Ik' spalloc(N*p,(N+1)*p,0)];
    end
    d2c = sparse(full(d2c));
elseif P==11
    d2c = zeros(4,4);
    H1_11 = -0.9*x(2)*0.25*exp(-0.5*x(1));
    H1_12 = -0.9*0.5*(1-exp(-0.5*x(1)));
    H2_11 = -0.8*x(2)*0.16*exp(-0.4*x(1));
    H2_12 = -0.8*0.4*(1-exp(-0.4*x(1)));
    d2c(1,1) = lam(1)*H1_11 - lam(2)*H1_11 + lam(3)*H2_11 - lam(4)*H2_11;
    d2c(1,2) = lam(1)*H1_12 - lam(2)*H1_12 + lam(3)*H2_12 - lam(4)*H2_12;
    d2c(2,1) = d2c(1,2);
elseif P==12
    d2c = zeros(5,5);
    d2c(1,1) = 2*sum(lam(1:3));
    d2c(2,2) = 2*sum(lam(1:3));
elseif P==13
    d2c = zeros(11,11);
    H9_9  =  (1-x(9))^(-2);
    H10_10 = (1-x(10))^(-2);
    H11_11 = (1-x(11))^(-2);
    d2c(9,9)   = lam(1)*H9_9;
    d2c(10,10) = lam(2)*H10_10;
    d2c(11,11) = lam(3)*H11_11;
elseif P==14
    d2c = zeros(5,5);
    d2c(4,4) = lam(1)*0.5*x(4)^(-1.5)*x(5)^2;
    d2c(4,5) = lam(1)*(-x(4)^(-0.5)*2*x(5));
    d2c(5,4) = d2c(4,5);
    d2c(5,5) = lam(1)*(4+0.5*x(5)^(-1.5)-4*x(4)^0.5);
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 080603 med  Switched to minlpAssign, cleaned