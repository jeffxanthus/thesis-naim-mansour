% geno_f.m
%
% function f = geno_f(x, Prob)
%
% Evaluates the objective function for test functions from the GENO
% development manual.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function f = geno_f(x, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        f = 5.3578547*x(3)^2 + 0.8356891*x(1)*x(5) + 37.293239*x(1) - 40792.141;
    case 2 % GENO The Colville #4 Function
        f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2 + 90*(x(4) - x(3)^2)^2 ...
            + (1 - x(3))^2 + 10.1*((x(2) - 1)^2 + (x(4) - 1)^2) + 19.8*(x(2) - 1)*(x(4) - 1);
    case 3 % GENO The 2-Dimensional Rastrigin Function
        f = x(1)^2 + x(2)^2 - cos(18*x(1)) - cos(18*x(2));
    case 4 % GENO The Economic Dispatch Problem
        f = 3*x(1) + 0.000001*x(1)^3 + 2*x(2) + 0.000002/(3*x(2)^3);
    case 5 % GENO A Pressure Vessel Design Problem
        x(1) = 0.0625*x(1);
        x(2) = 0.0625*x(2);
        f = 0.6224*x(1)*x(3)*x(4) + 1.7781*x(2)*(x(3)^2) + 3.1661*x(4)*(x(1)^2) ...
            + 19.84*x(3)*(x(1)^2);
    case 6 % GENO The Alkylation Process
        f = 5.04*x(1) + 0.035*x(2) + 10*x(3) + 3.36*x(5) - 0.063*x(4)*x(7);
    case 7 % GENO Decentralised Economic Planning
        f = (x(1)-10)^2 + 5*(x(2)-12)^2 + x(3)^4 + 3*(x(4)-11)^2 ...
            + 10*x(5)^6 + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);
    case 8 % GENO Heat Exchanger Optimisation
        f = x(1)+x(2)+x(3);
    case 9 % GENO The Harvest Problem
        f = -sum(sqrt(x(1:4)));
    case 10 % GENO A Non-linear Resource Allocation Problem
        f = 1;
        for i=1:5
            f = f*(1+i*x(i));
        end
        f = -f;
    case 11 % GENO Oligopolist Market Equilibrium Problem
        ci = [10;8;6;4;2];
        Ki = [5;5;5;5;5];
        alphai = [1/1.2; 1/1.1; 1.00; 1/0.9; 1/0.8];
        beta = 1/1.1;
        Q = sum(x);
        pQ = (5000/Q).^(beta);
        fq = ci.*x+(x.^(1+alphai)./(1+alphai)./Ki.^alphai);
        f = -(x.*pQ-fq);  % Maximise
    case 12 % GENO The Euclidean Compromise Solution I
        f1 = x(1)^2;
        f2 = (x(1)-2)^2;
        f = [f1;f2];
    case 13 % GENO The Euclidean Compromise Solution II
        f1 = (x(1)-1)^2+(x(2)-3)^2;
        f2 = (x(1)-4)^2+(x(2)-2)^2;
        f = [f1;f2]; % 0.5*sum(f1^2+f2^2);
    case 14 % GENO Efficient Portfolio Selection
        corr = [1 0.944 0.146 0.231 0.379 0.258; ...
            0.944 1 0.109 0.239 0.413 0.223; ...
            0.146 0.109 1 -0.169 -0.229 0.691; ...
            0.231 0.239 -0.169 1 0.882 -0.256; ...
            0.379 0.413 -0.229 0.882 1 -0.284; ...
            0.258 0.223 0.691 -0.256 -0.284 1];
        meanret = [15.852; 14.262; 31.336; 25.775; 50.228; 14.842];
        stddev  = [37.215; 41.773; 33.165; 62.009; 60.720; 23.757];
        f1 = sum(meanret.*x(1:6));
        f22 = 0;
        T = 6;
        f21 = x(1)^2*stddev(1)^2;
        for i=2:T
            f21 = f21 + x(i)^2*stddev(i)^2;
            for j=1:i-1
                f22 = f22 + corr(i,j)*x(i)*x(j)*stddev(i)*stddev(j);
            end
        end
        f2 = f21+2*f22;
        f2 = sqrt(f2);
        f = -f1+f2;
    case 15 % GENO A Dynamic Non-cooperative Game 1
        % u1(1) ... u1(T) u2(1) ... u2(T) .. x1(1) x1(T+1) x2(1) x2(T+1)
        T = 5;
        N = 2;
        f = x(T*N+T+1)-1/(2*T)*sum(x(T+1:2*T).^2);
        f = -f;
    case 16 % GENO A Dynamic Non-cooperative Game 2
        T = 5;
        N = 2;
        f = abs(x(T*N+T+1)) + abs(x(end));
        f = f + sum( abs(x(N*T+1:N*T+T)) + abs(x(N*T+2+T:end-1)) + 0.5*abs(x(1:T)) + 0.5*abs(x(T+1:T*N)) );
    case 17 % GENO A Dynamic Non-cooperative Game 3
        s = 0.5;
        q = 0.95;
        T = 5;
        v = 0.5;
        f = s*(x(end)^(1-v))*(q^T);
        f = f + sum( q.^(0:T-1)'.*x(1:T).^(1-v) );
        f = -f;
    case 18 % GENO A Dynamic Non-cooperative Game 4
        q = 1;
        s = 1;
        r = 1;
        T = 5;
        f = q*(x(end)^2);
        f = f + sum( s*x(T+1:end-1).^2 + r*x(1:T).^2);
    case 19 % GENO A Dynamic Non-cooperative Game 5
        f = 1;
        T = 5;
        for i=1:T
            f = f*(1+i*x(i));
        end
        f = -f;
    case 20 % GENO A Dynamic Non-cooperative Game 6
        T = 5;
        f = abs(x(end)-1) + sum( (2+x(1:T)).*3.^(-x(T+1:end-1)) );
    case 21 % GENO A Dynamic Non-cooperative Game 7
        T = 5;
        f = -sum(sqrt(x(1:T)));
end

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 080603  med  Switched to *Assign, cleaned