% function J = ls_J(x, Prob)
%
% Computes the Jacobian to least squares problem in the point x
% for the test problem P (Prob.P).
%
% J(i,j) is dr_i/d_x_j

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function J = ls_J(x, Prob)

global LS_A

x=x(:);
P=Prob.P;
uP=Prob.uP;
y=Prob.LS.y;
t=Prob.LS.t;
m=size(y,1);

if P==1
    % Powell
    J=[1; 2*uP(1)*x(1)+1];
elseif P==2
    % Walsh
    if x(1)==0
        b=1E100;
    else
        b=1/(x(1)*uP(1))-1;
    end
    J=zeros(m,2*size(x,2));
    for i=1:m
        a=1-x(1)*t(i)/x(2);
        if a <= 0, loga=1E100; else loga=log(a);end
        if x(2)==0
            J(i,1)=1E100;
            J(i,2)=1E100;
        elseif x(1)==0
            J(i,1)=1E100;
            J(i,2)=a^(b-1)*(x(1)*t(i)/x(2)^2)*b;
        else
            J(i,1)=a^b*(-1/(uP(1)*x(1)^2)*loga-b*t(i)/(x(2)*a));
            J(i,2)=a^(b-1)*(x(1)*t(i)/x(2)^2)*b;
        end
    end
elseif P==3
    % Gisela: uP(1) = K
    a=uP(1)*x(1)/(x(3)*(x(1)-x(2)));
    b=x(1)-x(2);
    J=zeros(m,3*size(x,2));
    for i=1:m
        e1=exp(-x(1)*t(i)); e2=exp(-x(2)*t(i));
        J(i,1)=a*(t(i)*e1+(e2-e1)*(1-1/b));
        J(i,2)=a*(-t(i)*e2+(e2-e1)/b);
        J(i,3)=-a*(e2-e1)/x(3);
    end
elseif P==4
    m=round(length(y(:))/2);
    J=zeros(2*m,3*size(x,2));
    if isempty(LS_A)
        r = nlp_r(x,Prob);
    end
    J(1:2*m,1:3)=[LS_A,zeros(2*m,2)];
    J(1:m,2)=ones(m,1);
    J(m+1:2*m,3)=ones(m,1);
elseif P==5
    % Population problem
    J=zeros(m,2*size(x,2));
    for i=1:m, J(i,:)=[x(2)^t(i),t(i)*x(1)*x(2)^(t(i)-1)]; end
elseif P==6
    % Plasmid n=2
    p_0 = uP(2);
    J=zeros(m,2*size(x,2));
    for i=1:m
        e1 = exp((x(1)+x(2))*t(i));
        p1 = (p_0*x(1)+x(2));
        j11 = p_0*e1+p1*t(i)*e1;
        j12 = p1*e1+x(1)*(1-p_0);
        j13 = (p1*e1-x(2)*(1-p_0));
        j14 = (p_0*e1+p1*t(i)*e1+1-p_0);
        j15 = (p1*e1+x(1)*(1-p_0))^2;
        j21 = e1+p1*t(i)*e1-1+p_0;
        j24 = (e1+p1*t(i)*e1);
        J(i,1) = j11/j12 - j13*j14/j15;
        J(i,2) = j21/j12 - j13*j24/j15;
    end
elseif P==7
    % Plasmid n=3
    J=zeros(m,3*size(x,2));
    for i=1:m
        e1 = exp((x(1)+x(2))*t(i));
        p1 = (x(3)*x(1)+x(2));
        j11 = x(3)*e1+p1*t(i)*e1;
        j12 = p1*e1+x(1)*(1-x(3));
        j13 = (p1*e1-x(2)*(1-x(3)))*(x(3)*e1+p1*t(i)*e1+1-x(3));
        j14 = (p1*e1+x(1)*(1-x(3)))^2;
        j21 = e1+p1*t(i)*e1-1+x(3);
        j22 = p1*e1+x(1)*(1-x(3));
        j23 = (p1*e1-x(2)*(1-x(3)))*(e1+p1*t(i)*e1);
        j24 = (p1*e1+x(1)*(1-x(3)))^2;
        j31 = x(1)*e1+x(2);
        j32 = p1*e1+x(1)*(1-x(3));
        j33 = (p1*e1-x(2)*(1-x(3)))*(x(1)*e1-x(1));
        j34 = (p1*e1+x(1)*(1-x(3)))^2;
        J(i,1) = j11/j12-j13/j14;
        J(i,2) = j21/j22-j23/j24;
        J(i,3) = j31/j32-j33/j34;
    end
elseif P==8
    % Plasmid n=3 (subst.)
    D   = uP(3);
    J=zeros(m,3*size(x,2));
    for i=1:m
        e1 = exp((x(1)-x(2)+x(3))*t(i));
        d1 = (D-x(2)+x(3));
        j11 = x(3)*(D-x(2))/(x(1)-x(2))^2+1;
        j12 = d1*e1+x(1)-D;
        j13 = (x(3)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*(d1*t(i)*e1+1);
        j14 = (d1*e1+x(1)-D)^2;
        j21 = x(3)*(1/(x(1)-x(2))-(D-x(2))/(x(1)-x(2))^2);
        j22 = d1*e1+x(1)-D;
        j23 = (x(3)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*(-e1-d1*t(i)*e1);
        j24 = (d1*e1+x(1)-D)^2;
        j31 = 1-(D-x(2))/(x(1)-x(2));
        j32 = d1*e1+x(1)-D;
        j33 = (x(3)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*(e1+d1*t(i)*e1);
        j34 = (d1*e1+x(1)-D)^2;
        J(i,1) = -j11/j12 + j13/j14;
        J(i,2) = -j21/j22 + j23/j24;
        J(i,3) = -j31/j32 + j33/j34;
    end
elseif P==9
    % Plasmid n=3 (probability)
    D   = uP(3);
    J=zeros(m,3*size(x,2));
    for i=1:m
        e1 = exp((x(1)+x(2)*(x(3)-1))*t(i));
        d1 = (D+x(2)*(x(3)-1));
        j11 = x(3)*x(2)*(D-x(2))/(x(1)-x(2))^2+1;
        j12 = d1*e1+x(1)-D;
        j13 = (x(3)*x(2)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*(d1*t(i)*e1+1);
        j14 = (d1*e1+x(1)-D)^2;
        j21 = x(3)*(1-(D-x(2))/(x(1)-x(2)))+x(3)*x(2)*(1/(x(1)-x(2))-(D-x(2))/(x(1)-x(2))^2);
        j22 = d1*e1+x(1)-D;
        j23 = (x(3)*x(2)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*((x(3)-1)*e1+d1*(x(3)-1)*t(i)*e1);
        j24 = (d1*e1+x(1)-D)^2;
        j31 = x(2)*(1-(D-x(2))/(x(1)-x(2)));
        j32 = d1*e1+x(1)-D;
        j33 = (x(3)*x(2)*(1-(D-x(2))/(x(1)-x(2)))+x(1)-D)*(x(2)*e1+d1*x(2)*t(i)*e1);
        j34 = (d1*e1+x(1)-D)^2;
        J(i,1) = -j11/j12+j13/j14;
        J(i,2) = -j21/j22+j23/j24;
        J(i,3) = -j31/j32+j33/j34;
    end
elseif P==10
    % Parameterized test function (Huschens)
    phi = uP(1);
    J=zeros(m,2*size(x,2));
    J(1,1) = 1;
    J(2,1) = x(2);
    J(2,2) = x(1)-2*phi;
    J(3,2) = 1;
elseif P==11
    % Signomial problem
    n = uP(1);
    l = Prob.user.l;
    C = Prob.user.C;
    A = Prob.user.A;
    J = zeros(m,n);
    for i=1:m
        for j=1:n
            for k=1:l
                z=C(i,k);
                for t=1:n
                    if t~=j
                        z=z*x(t)^A(i,t,k);
                    else
                        z=z*A(i,j,k)*x(j)^(A(i,j,k)-1);
                    end
                end
                J(i,j)=J(i,j)+z;
            end
        end
    end

elseif P==12
    % Signomial problem
    n = uP(1);
    l = Prob.user.l;
    C = Prob.user.C;
    A = Prob.user.A;
    J = zeros(m,n);
    for i=1:m
        for j=1:n
            for k=1:l
                z=C(i,k);
                for t=1:n
                    if t~=j
                        z=z*x(t)^A(i,t,k);
                    else
                        z=z*A(i,j,k)*x(j)^(A(i,j,k)-1);
                    end
                end
                J(i,j)=J(i,j)+z;
            end
        end
    end
elseif P==13
    % Exponential problem (rand)
    n = uP(1);
    l = Prob.user.l;
    C = Prob.user.C;
    A = Prob.user.A;
    J = zeros(m,n);
    for i=1:m
        for j=1:n
            for k=1:l
                z=0;
                for t=1:n
                    z=z+x(t)*A(i,t,k);
                end
                J(i,j)=J(i,j)+C(i,k)*exp(z)*A(i,j,k);
            end
        end
    end
elseif P==14
    % Exponential problem (pseudorand)
    n = uP(1);
    l = Prob.user.l;
    C = Prob.user.C;
    A = Prob.user.A;
    J = zeros(m,n);
    for i=1:m
        for j=1:n
            for k=1:l
                z=0;
                for t=1:n
                    z=z+x(t)*A(i,t,k);
                end
                J(i,j)=J(i,j)+C(i,k)*exp(z)*A(i,j,k);
            end
        end
    end
elseif P==15
    % Trigonometric problem
    n = uP(1);
    A = Prob.user.A;
    B = Prob.user.B;
    C = Prob.user.C;
    J = zeros(m,n);
    r =-C+A*sin(x)+B*cos(x);
    for i=1:m
        for j=1:n
            J(i,j)=2*r(i)*(A(i,j)*cos(x(j))-B(i,j)*sin(x(j)));
        end
    end
end

% MODIFICATION LOG:
%
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981108  hkh  Use general global variable LS_A instead of circle_t
% 981210  hkh  Safeguard Walsh against division by 0 and log of negative numb
% 030323  hkh  Adding four problems
% 031201  hkh  Revising AD handling, new for MAD, changes for ADMAT
% 031201  hkh  Use zeros(m,n*size(x,2)) to make AD tool MAD work, size(x,2)=1
% 040414  hkh  Check if LS_A empty, and call nlp_r if so
% 050302  hkh  Added problems 11-15
% 080603  med  Switched to clsAssign, cleaned