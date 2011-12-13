% mgh_r.m
%
% function r = mgh_r(x, Prob)
%
% Computes residuals to More, Garbow, Hillstrom nonlinear least
% squares problem in the point x for the test problem Prob.P

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function r = mgh_r(x, Prob)

P=Prob.P;
y=Prob.LS.y(:);
t=Prob.LS.t;

m=size(y,1);
n=Prob.N;

if P==1
    % Rosenbrock [More G H #1]
    r = zeros(m,1);
    r(1) = 10 * (x(2) - x(1)^2);
    r(2) = 1 - x(1);
elseif P==2
    % Freudenstein and Roth [More G H #2]
    r = zeros(m,1);
    r(1) = -13 + x(1) + ((5 - x(2))*x(2) - 2)*x(2);
    r(2) = -29 + x(1) + (( x(2) + 1)*x(2) - 14)*x(2);
elseif P==3
    % Powell Badly Scaled [More G H #3]
    r = zeros(m,1);
    r(1) =  1e4*x(1)*x(2) - 1;
    r(2) =  exp(-x(1)) + exp(-x(2)) - 1.0001;
elseif P==4
    % Brown Badly Scaled [More G H #4]
    r = zeros(m,1);
    r(1) =  x(1) - 1e6;
    r(2) =  x(2) - 2e-6;
    r(3) =  x(1)*x(2) - 2;
elseif P==5
    % Beale [More G H #5]
    r = -x(1)*(1 - x(2).^(1:3)');
elseif P==6
    % Jennrich and Sampson [More G H #6]
    r = 2 + 2*(1:m)' - exp(x(1)*(1:m)') - exp(x(2)*(1:m)');
elseif P==7
    % Helical Valley [More G H #7]
    r = zeros(m,1);
    if x(1)==0 & x(2)~=0
        quote=Inf;
    elseif x(1)==0 & x(2)==0
        quote=1;
    else
        quote=x(2)/x(1);
    end
    r(1) = 10*( x(3) - 10*(0.5*atan( quote ) / pi + 0.5*(x(1)<0)) );
    r(2) = 10*(norm(x(1:2))-1);
    r(3) = x(3);
elseif P==8
    % Bard [More G H #8]
    if x(2)==0 & x(3)==0
        r=Inf*ones(m,1);
    else
        r = x(1)+ (1:m)'./((16-(1:m)')*x(2) + min((1:m),16-(1:m))'*x(3));
    end
elseif P==9
    % Gaussian [More G H #9]
    r = x(1).*exp((-x(2).*(t-x(3)).^2)./2);
elseif P==10
    % Meyer [More G H #10]
    r = x(1)*exp( x(2)./(t + x(3)) );
elseif P==11
    % Gulf Research and Development [More G H #11]
    if x(1)==0
        quote=Inf;
    else
        quote=1/x(1);
    end
    r = exp( -abs( 25 + (-50.*log(t) ).^ (2/3) - x(2)).^x(3)*quote ) - t;
elseif P==12
    % Box 3-Dimensional [More G H #12]
    r = exp(-t.*x(1)) - exp(-t.*x(2)) - x(3)*( exp(-t) - exp(-t.*10) );
elseif P==13
    % Powell Singular [More G H #13]
    r = zeros(m,1);
    r(1) = x(1) + 10*x(2);
    r(2) = sqrt(5) * (x(3) - x(4));
    r(3) = (x(2) - 2*x(3))^2;
    r(4) = sqrt(10) * (x(1) - x(4))^2;
elseif P==14
    % Wood [More G H #14]
    r = zeros(m,1);
    r(1) = 10*(x(2) - x(1)*x(1));
    r(2) = 1 - x(1);
    r(3) = sqrt(90)*(x(4) - x(3)*x(3));
    r(4) = 1 - x(3);
    r(5) = sqrt(10)*(x(2) + x(4) - 2);
    r(6) = (x(2) - x(4))/sqrt(10);
elseif P==15
    % Kowalik and Osborne [More G H #15]
    u = [4 2 1 0.5 0.25 0.167 0.125 0.1 0.0833 0.0714 0.0625]';
    r = x(1).*u.*(u + x(2))./(u.*(u + x(3)) + x(4));
elseif P==16
    % Brown and Dennis (More,Garbow,Hillstrom #16)
    r = (x(1) + t*x(2) - exp(t)).^2 + (x(3) + x(4)*sin(t) - cos(t)).^2;
elseif P==17
    % Osborne1 (More,Garbow,Hillstrom #17)
    r = x(1) + x(2).*exp(-t.*x(4)) + x(3).*exp(-t.*x(5));
elseif P==18
    % Biggs (More,Garbow,Hillstrom #18)
    r = x(3)*exp(-t.*x(1)) + x(6)*exp(-t.*x(5)) - x(4)*exp(-t.*x(2));
elseif P==19
    % Osborne2 (More,Garbow,Hillstrom #19)
    r = x(1)*exp(-t.*x(5)) + x(2)*exp(-(t - x(9)).^2*x(6)) ...
        + x(3)*exp(-(t - x(10)).^2*x(7)) + x(4)*exp(-(t - x(11)).^2*x(8));
elseif P==20
    % Watson [More G H #20]'
    r = zeros(m,1);
    for i = 1:m-2
        s2 = (t(i).^(0:n-1))*x;
        s1 = (1:n-1)*(  x(2:n).*(t(i).^(0:n-2))' );
        r(i) =  s1 - s2*s2 - 1;
    end
    r(30) = x(1);
    r(31) = x(2) - x(1)*x(1) - 1;
elseif P==21
    % Extended Rosenbrock [More G H #21]
    r = zeros(m,1);
    r(1:2:n-1) = 10 * (x(2:2:n) - x(1:2:n-1).^2);
    r(2:2:n)   = 1 - x(1:2:n-1);
elseif P==22
    % Extended Powell Singular [More G H #22]
    r = zeros(m,1);
    r(1:4:n) = x(1:4:n) + 10*x(2:4:n);
    r(2:4:n) = sqrt(5) * (x(3:4:n) - x(4:4:n));
    r(3:4:n) = ( x(2:4:n) - 2*x(3:4:n) ).^2;
    r(4:4:n) = sqrt(10)*(x(1:4:n) - x(4:4:n)).^2;
elseif P==23
    % Penalty I [More G H #23]
    r = zeros(m,1);
    r(1:n) = sqrt(1e-5)*(x - 1);
    r(m) = sum(x.^2) - 0.25;
elseif P==24
    % Penalty II [More G H #24]
    r = zeros(m,1);
    r(1)    = x(1) - 0.2;
    r(2:n) = sqrt(1e-5)*( ( exp(x(2:n)*0.1) + exp(x(1:n-1)*0.1)) ...
        - (exp((2:n)'*0.1) + exp((1:n-1)'*0.1)) );
    r(n+1:2*n-1) = sqrt(1e-5)*( exp(x(2:n)*0.1)  - exp(-0.1) );
    r(m) = (n-(1:n)+1)*x.^2 -1;
elseif P==25
    % Variably Dimensioned [More G H #25]
    r = zeros(m,1);
    summ = (1:n)*(x-1);
    r(1:n) = x-1;
    r(n+1) = summ;
    r(n+2) = summ^2;
elseif P==26
    % Trigonometric [More G H #26]
    r = n + (1:n)'.*(1 - cos(x)) - sin(x) - sum(cos(x));
elseif P==27
    % Brown Almost Linear [More G H #27]
    r = x + sum(x) - (n+1) ;
    r(n) = prod(x)-1;
elseif P==28
    % Discrete Boundary Value [More G H #28]
    xm1 = [0;x(1:n-1)];
    xp1 = [x(2:n);0];
    r = 2*x - xm1 - xp1 + (x+t+1).^3/( 2*(n+1)^2 );
elseif P==29
    % Discrete Integral Equation [More G H #29]
    r = x;
    h  = 1/(n+1);
    h2 = 0.5*h;
    for j = 1:n
        xj = x(j);
        tj = j*h;
        oj = 1 - tj;
        bj = (xj + tj + 1);
        sj = h2*bj*bj;
        cj = sj*bj;
        for i = 1:n
            ti = i*h;
            oi = 1 - ti;
            if (j <= i), r(i) = r(i) + oi*tj*cj; end
            if (j > i),  r(i) = r(i) + ti*oj*cj; end
        end
    end
elseif P==30
    % Broyden Tridiagonal [More G H #30]
    xm1 = [0;x(1:n-1)];
    xp1 = [x(2:n);0];
    r = (3-2*x).*x - xm1 - 2*xp1 + 1;
elseif P==31
    % Broyden Banded [More G H #31]
    r = zeros(m,1);
    ilbw   = 5;
    iubw   = 1;
    for i = 1:n
        xi      = x(i);
        t       = 5*xi*xi;
        r(i)    = xi*(2 + t) + 1;
        j1      = max(1,i-ilbw);
        j2      = min(n,i+iubw);
        for j = j1:j2
            xj = x(j);
            if (j ~= i),  r(i)   = r(i) - xj*(1 + xj); end
        end
    end
elseif P==32
    % Linear --- Full Rank [More G H #32]
    r = [x-2*sum(x)/m-1;ones(m-n,1)*(-2*sum(x)/m-1)];
elseif P==33
    % Linear --- Rank 1 [More G H #33]
    r = (1:m)'*( (1:n)*x ) -1;
elseif P==34
    % Linear --- Rank 1 with Zero Columns and Rows [More G H #34]
    r = zeros(m,1);
    r(1) = -1;
    r(2:m-1) = ((2:m-1)'-1)*( (2:n-1)*x(2:n-1) )-1;
    r(m) = -1;
elseif P==35
    % Chebyquad [More G H #35]
    r = zeros(m,1);
    coeff = 1/n;
    for i = 1:m/2
        i2 = 2*i;
        r(i2) = 1/((i2*i2) - 1);
    end
    for j = 1:n
        t0 = 1;
        s0 = 0;
        s1 = 2;
        t1 = 2 * x(j) - 1;
        tj = t1;
        r(1) = r(1) + coeff*t1;
        for i = 2:m
            u  = 2*tj;
            t  = u*t1 - t0;
            s  = 4*t1 + u*s1 - s0;
            r(i) = r(i) + coeff*t;
            t0 = t1;
            t1 = t;
            s0 = s1;
            s1 = s;
        end
    end
end

if Prob.LS.yUse & m==length(r),  r=r-y; end

% MODIFICATION LOG:
%
% 981018  hkh  Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981118  hkh  Use new flag Prob.NLLS.UseYt
% 990314  hkh  Remove unnecessary yMod computation
% 050303  hkh  Use n=Prob.N instead of uP
% 080603  med  Switched to clsAssign, cleaned