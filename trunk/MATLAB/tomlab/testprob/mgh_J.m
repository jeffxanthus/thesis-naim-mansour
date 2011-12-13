% mgh_J.m
%
% function J = mgh_J(x, Prob)
%
% Computes the Jacobian to More, Garbow, Hillstrom nonlinear least
% squares test problems in the point x for the test problem Prob.P
%
% J(i,j) is dr_i/d_x_j

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function J = mgh_J(x, Prob)

x=x(:);

P=Prob.P;
y=Prob.LS.y;
t=Prob.LS.t;
m=size(y,1);
n=Prob.N;

if P==1
    % Rosenbrock [More G H #1]
    J = [-20 * x(1),10;-1,0];
elseif P==2
    % Freudenstein and Roth [More G H #2]
    J=zeros(m,2);
    J(1,1) = 1;
    J(2,1) = 1;
    J(1,2) = (10 - 3*x(2))*x(2) - 2;
    J(2,2) = (3*x(2) + 2)*x(2) - 14;
elseif P==3
    % Powell Badly Scaled [More G H #3]
    J=zeros(m,2);
    J(1,1) = 1e4*x(2);
    J(2,1) = -exp(-x(1));
    J(1,2) = 1e4*x(1);
    J(2,2) = -exp(-x(2));
elseif P==4
    % Brown Badly Scaled [More G H #4]
    J = [1 0;0 1;x(2) x(1)];
elseif P==5
    % Beale [More G H #5]
    J=zeros(m,2);
    for i = 1:3
        p = x(2)^(i-1);
        s = p*x(2);
        J(i,1) = s - 1;
        J(i,2) = x(1)*i*p;
    end
elseif P==6
    % Jennrich and Sampson [More G H #6]
    J=zeros(m,2);
    for i = 1:m
        J(i,1) = -i*exp(i*x(1));
        J(i,2) = -i*exp(i*x(2));
    end
elseif P==7
    % Helical Valley [More G H #7]
    J=zeros(m,3);
    piinv  = 0.25 / atan(1);
    t = x(1)*x(1) + x(2)*x(2);
    if t==0
        q=Inf;
    else
        q = 1 / t;
    end
    p = sqrt(q);
    fpq = 50*piinv*q;
    tr  = 10*p;
    J(1,1) =  fpq*x(2);
    J(1,2) = -fpq*x(1);
    J(1,3) =  10;
    J(2,1) =  tr*x(1);
    J(2,2) =  tr*x(2);
    J(2,3) =  0;
    J(3,1) =  0;
    J(3,2) =  0;
    J(3,3) =  1;
elseif P==8
    % Bard [More G H #8]
    J=zeros(m,3);
    for i = 1:m
        ui = i;
        vi = 16-i;
        wi = min(ui,vi);
        di = vi*x(2) + wi*x(3);
        if di==0
            qi = Inf;
        else
            qi = -ui / (di*di);
        end
        J(i,1) = 1;
        J(i,2) = qi * vi;
        J(i,3) = qi * wi;
    end
elseif P==9
    % Gaussian [More G H #9]
    J=zeros(m,3);
    for i = 1:m
        t3 = t(i) - x(3);
        tt = (t3*t3)/2;
        ei = exp(-x(2)*tt);
        si = x(1)*ei;
        J(i,1) =  ei;
        J(i,2) = -tt*si;
        J(i,3) =  x(2)*t3*si;
    end
elseif P==10
    % Meyer [More G H #10]
    d1 = 1./(t+x(3));
    e1 =  exp(x(2)*d1);
    J = [e1, x(1)*d1.*e1,-x(1)*x(2)*d1.^2.*e1];
elseif P==11
    % Gulf Research and Development [More G H #11]
    J=zeros(m,3);
    if x(1)==0
        quote=Inf;
    else
        quote=1/x(1);
    end
    for i = 1:m
        yi = 25 + (-50*log(t(i))) ^ (2/3);
        ai = yi - x(2);
        av = abs(ai);
        bi = av^x(3);
        ci = bi*quote;
        ei = exp(-ci);
        d1 =   ci * quote;
        d2 =   x(3) * quote * av^(x(3) - 1);
        d3 = - log(av) * ci;
        J(i,1) = d1 * ei;
        if (ai >= 0), J(i,2) =  d2 * ei; end
        if (ai <  0), J(i,2) = -d2 * ei; end
        J(i,3) = d3 * ei;
    end
elseif P==12
    % Box 3-Dimensional [More G H #12]
    J=zeros(m,3);
    for i = 1:10
        e1 = exp(-t(i)*x(1));
        e2 = exp(-t(i)*x(2));
        t3 = exp(-t(i)) - exp(-t(i)*10);
        J(i,1) = -t(i)*e1;
        J(i,2) =  t(i)*e2;
        J(i,3) = -t3;
    end
elseif P==13
    % Powell Singular [More G H #13]
    J=zeros(m,4);
    t1 = x(2) - 2*x(3);
    t0 = x(1) - x(4);
    s0 = 2*sqrt(10) * t0;
    J(1,1) =  1;
    J(1,2) =  10;
    J(2,3) =  sqrt(5);
    J(2,4) = -sqrt(5);
    J(3,2) =  2 * t1;
    J(3,3) = -4 * t1;
    J(4,1) =  s0;
    J(4,4) = -s0;
elseif P==14
    % Wood [More G H #14]
    J = [-20*x(1),10,0,0;-1,0,0,0;0,0,-2*sqrt(90)*x(3),sqrt(90);0,0,-1,0 ...
        ;0,sqrt(10),0,sqrt(10);0,1/sqrt(10),0,-1/sqrt(10) ];
elseif P==15
    % Kowalik and Osborne [More G H #15]
    J=zeros(m,4);
    u = [4 2 1 0.5 0.25 0.167 0.125 0.1 0.0833 0.0714 0.0625]';
    for i = 1:m
        t2 = u(i)*(u(i) + x(2));
        di = u(i)*(u(i) + x(3)) + x(4);
        qi = 1 / di;
        p  =  t2*qi;
        s  =  x(1)*qi;
        v  = -p*s;
        J(i,1) = p;
        J(i,2) = s*u(i);
        J(i,3) = v*u(i);
        J(i,4) = v;
    end
elseif P==16
    % Brown and Dennis (More,Garbow,Hillstrom #16)
    J = [2*(x(1) + t.*x(2) - exp(t)) , (2*(x(1) + t.*x(2) - exp(t))).*t ...
        2*(x(3) + x(4).*sin(t) - cos(t)) , 2*(x(3) + x(4).*sin(t) - cos(t)).*sin(t)];
elseif P==17
    % Osborne1 (More,Garbow,Hillstrom #17)
    J=zeros(m,5);
    for i = 1 : 33
        e4 = exp(-t(i)*x(4));
        e5 = exp(-t(i)*x(5));
        t2 = x(2)*e4;
        t3 = x(3)*e5;
        J(i,1) =  1;
        J(i,2) =  e4;
        J(i,3) =  e5;
        J(i,4) = -t(i)*t2;
        J(i,5) = -t(i)*t3;
    end
elseif P==18
    % Biggs (More,Garbow,Hillstrom #18)
    e1 = exp(-t*x(1));
    e2 = exp(-t*x(2));
    e5 = exp(-t*x(5));
    J = [-x(3)*t.*e1, x(4)*t.*e2, e1,-e2, -x(6)*t.*e5, e5];
elseif P==19
    % Osborne2 (More,Garbow,Hillstrom #19)
    J=zeros(m,11);
    for i = 1 : 65
        t9   = t(i) - x(9);
        t10  = t(i) - x(10);
        t11  = t(i) - x(11);
        s9   = t9 * t9;
        s10  = t10 * t10;
        s11  = t11 * t11;
        e1   = exp(-t(i)*x(5));
        e2   = exp(-s9*x(6));
        e3   = exp(-s10*x(7));
        e4   = exp(-s11*x(8));
        r2   = x(2)*e2;
        r3   = x(3)*e3;
        r4   = x(4)*e4;
        J(i,1) = e1;
        J(i,2) = e2;
        J(i,3) = e3;
        J(i,4) = e4;
        J(i,5) = -t(i)*x(1)*e1;
        J(i,6) = -s9*r2;
        J(i,7) = -s10*r3;
        J(i,8) = -s11*r4;
        J(i,9) =  2*t9*x(6)*r2;
        J(i,10) =  2*t10*x(7)*r3;
        J(i,11) =  2*t11*x(8)*r4;
    end
elseif P==20
    % Watson [More G H #20]
    J=zeros(m,n);
    n1 = n - 1;
    n2 = n1 - 1;
    r2 = x(n);
    r1 = n1*r2;
    for i = 1:m-2
        s1 = r1;
        s2 = r2;
        nj = n1;
        nj1= n2;
        for j = 2:n1
            xnj = x(nj);
            s1  = nj1*xnj + t(i)*s1;
            s2  = xnj + t(i)*s2;
            nj  = nj1;
            nj1 = nj1 - 1;
        end
        s2   =  x(1) + t(i)*s2;
        t2   =  2*s2;
        J(i,1)  = -t2;
        t2i     =  t2*t(i);
        J(i,2)  =  1 - t2i;
        jm1     = 2;
        tijm2   = t(i);
        for j = 3:n
            J(i,j) = (jm1 - t2i) * tijm2;
            jm1    = j;
            tijm2  = tijm2 * t(i);
        end
    end
    for j = 3:n
        J(30,j) = 0;
        J(31,j) = 0;
    end
    J(30,1) =  1;
    J(31,1) = -2*x(1);
    J(30,2) =  0;
    J(31,2) =  1;
elseif P==21
    % Extended Rosenbrock [More G H #21]
    J=zeros(m,n);
    for i = 1:n/2
        i2    = 2*i;
        i2m1  = i2 - 1;
        xi2m1 = x(i2m1);
        J(i2m1,i2m1) = -10 * 2 * xi2m1;
        J(i2m1,  i2) = 10;
        J(i2,i2m1)   = -1;
    end
elseif P==22
    % Extended Powell Singular [More G H #22]
    J=zeros(m,n);
    for i = 1:n/4;
        i4   = 4*i;
        i4m1 = i4 - 1;
        i4m2 = i4 - 2;
        i4m3 = i4 - 3;
        xi4   = x(i4);
        xi4m1 = x(i4m1);
        xi4m2 = x(i4m2);
        xi4m3 = x(i4m3);
        t1 = xi4m2 - 2*xi4m1;
        t0 = xi4m3 - xi4;
        s0 = 2*sqrt(10) * t0;
        J(i4m3,i4m3) =  1;
        J(i4m3,i4m2) =  10;
        J(i4m2,i4m1) =  sqrt(5);
        J(i4m2,i4)   = -sqrt(5);
        J(i4m1,i4m2) =  2 * t1;
        J(i4m1,i4m1) = -4 * t1;
        J(i4,i4m3)   =  s0;
        J(i4,i4)     = -s0;
    end
elseif P==23
    % Penalty I [More G H #23]
    J=zeros(m,n);
    a2 = sqrt(1e-5);
    summ = - 0.25;
    for j = 1:n
        xj   = x(j);
        summ  = summ + xj*xj;
        J(j,j) = a2;
        J(m,j) = 2*xj;
    end
elseif P==24
    % Penalty II [More G H #24]
    J=zeros(m,n);
    a     = 1.d-5;
    ar = sqrt(a);
    scale = 0.1;
    x1  = x(1);
    summ     = n*x1*x1;
    J(1,1)  = 1;
    J(m,1)  = 2*n*x1;
    tj  = ar*exp(x1*scale);
    uj  = scale*tj;
    jm1 = 1;
    for j = 2:n
        xj   = x(j);
        ujm1 = uj;
        tj   = ar*exp(xj*scale);
        uj   = scale*tj;
        vj   = (n-jm1)*xj;
        summ  = summ + vj*xj;
        J(j,jm1)   = ujm1;
        J(j,j)     = uj;
        J(n+jm1,j) = uj;
        J(m,j)     = 2*vj;
        jm1 = j;
    end
elseif P==25
    % Variably Dimensioned [More G H #25]
    J=zeros(m,n);
    summ = 0;
    for j = 1:n
        xjm1 = x(j) - 1;
        summ  = summ + j*xjm1;
    end
    for j = 1:n
        J(j,j)   = 1;
        J(n+1,j) = j;
        J(n+2,j) = 2*summ*j;
    end
elseif P==26
    % Trigonometric [More G H #26]
    J=zeros(m,n);
    for j = 1:n
        xj  = x(j);
        sxj = sin(xj);
        for i = 1:n
            J(i,j) = sxj;
        end
        J(j,j) = (j+1)*sxj - cos(xj);
    end
elseif P==27
    % Brown Almost Linear [More G H #27]
    J=zeros(m,n);
    for j = 1:n
        for i = 1:n-1
            J(i,j) = 1;
        end
        J(j,j) = 2;
    end
    for j = 1:n
        t = 1;
        for i = 1:n
            if (i ~= j)
                t = t * x(i);
            end
        end
        J(n,j) = t;
    end
elseif P==28
    % Discrete Boundary Value [More G H #28]
    J=zeros(m,n);
    i  = 1;
    xi = x(1);
    h  = 1/(n+1);
    im1=[];
    for ip1 = 2:n+1;
        if (i ~= n)
            xip1 = x(ip1);
        end
        ti = i*h;
        ri = (xi + ti + 1);
        si = h*h*ri*ri/2;
        J(i,i) = 2 + 3*si;
        if (i ~= 1)
            J(i,im1) = -1;
        end
        if (i ~= n)
            J(i,ip1) = -1;
        end
        im1  = i;
        i    = ip1;
        xi   =  xip1;
    end
elseif P==29
    % Discrete Integral Equation [More G H #29]
    h  = 1/(n+1);
    h2 = 0.5*h;
    J=zeros(m,n);
    for j = 1:n
        xj = x(j);
        tj = j*h;
        oj = 1 - tj;
        bj = (xj + tj + 1);
        sj = h2*bj*bj;
        sj = sj*3;
        J(j,j) = 1 + oj*tj*sj;
        for i = 1:n
            ti = i*h;
            oi = 1 - ti;
            if (j < i),  J(i,j) = oi*tj*sj; end
            if (j > i),  J(i,j) = ti*oj*sj; end
        end
    end
elseif P==30
    % Broyden Tridiagonal [More G H #30]
    J=zeros(m,n);
    i = 1;
    xi = x(1);
    for ip1 = 2:n+1
        if (i < n), xip1 = x(ip1); end
        J(i,i) = 3 - 4*xi;
        if (i > 1),  J(i,im1) = -1; end
        if (i < n),  J(i,ip1) = -2; end
        im1  = i;
        i    = ip1;
        xi   = xip1;
    end
elseif P==31
    % Broyden Banded [More G H #31]
    J=zeros(m,n);
    ilbw   = 5;
    iubw   = 1;
    for i = 1:n
        xi      = x(i);
        t       = 5*xi*xi;
        J(i,i)  = 2 + 3*t;
        j1      = max(1,i-ilbw);
        j2      = min(n,i+iubw);
        for j = j1:j2
            xj = x(j);
            if (j ~= i),  J(i,j) = -(1 + 2*xj); end
        end
    end
elseif P==32
    % Linear --- Full Rank [More G H #32]
    J=zeros(m,n);
    dnfi2 = 2/m;
    for i = 1:m
        for j = 1:n
            if (i <= n)
                J(i,j) = -dnfi2;
                J(i,i) = 1 - dnfi2;
            end
            if (i > n),  J(i,j) = -dnfi2; end
        end
    end
elseif P==33
    % Linear --- Rank 1 [More G H #33]
    J=zeros(m,n);
    summ = 0;
    for j = 1:n
        summ = summ + x(j)*j;
    end
    for i = 1:m
        for j = 1:n
            J(i,j) = i*j;
        end
    end
elseif P==34
    % Linear --- Rank 1 with Zero Columns and Rows [More G H #34]
    J=zeros(m,n);
    summ = 0;
    for j = 1:n
        if (j > 1) & (j < n), summ = summ + x(j)*j; end
    end
    for i = 2:m-1
        di = (i-1);
        for j = 2:n-1
            J(i,j) = di*j;
        end
    end
elseif P==35
    % Chebyquad [More G H #35]
    J=zeros(m,n);
    coeff = 1/n;
    for j = 1:n
        t0 = 1;
        s0 = 0;
        s1 = 2;
        t1 = 2 * x(j) - 1;
        tj = t1;
        J(1,j) = coeff*s1;
        for i = 2:m
            u  = 2*tj;
            t  = u*t1 - t0;
            s  = 4*t1 + u*s1 - s0;
            J(i,j) = coeff*s;
            t0 = t1;
            t1 = t;
            s0 = s1;
            s1 = s;
        end
    end
end

% MODIFICATION LOG
%
% 980911  mbk  f changed to J on the last line.
% 981011  hkh  Changed to use NLPLIB Gateway routines
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 990912  hkh  Correct error in J(3,3) in #7 Helical Valley.
% 990912  hkh  Change computation of #10 Meyer.
% 050303  hkh  Use n=Prob.N instead of uP
% 080603  med  Switched to clsAssign, cleaned