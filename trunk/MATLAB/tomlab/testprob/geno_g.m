% geno_g.m
%
% function g = geno_g(x, Prob)
%
% Evaluates the gradient for test functions from the GENO development
% manual.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function g = geno_g(x, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        g = zeros(5,1);
        g(1,1) = 37.293239 + 0.8356891*x(5);
        g(3,1) = 2*5.3578547*x(3);
        g(5,1) = 0.8356891*x(1);
    case 2 % GENO The Colville #4 Function
        g = zeros(4,1);
        g(1,1) = -2*100*(x(2) - x(1)^2)*2*x(1) - 2*(1 - x(1));
        g(2,1) = 2*100*(x(2) - x(1)^2) + 2*10.1*((x(2) - 1)) + 19.8*(x(4) - 1);
        g(3,1) = -2*90*(x(4) - x(3)^2)*2*x(3) - 2*(1 - x(3));
        g(4,1) = 2*90*(x(4) - x(3)^2) + 2*10.1*(x(4) - 1) + 19.8*(x(2) - 1);
    case 3 % GENO The 2-Dimensional Rastrigin Function
        g = zeros(2,1);
        g(1,1) = 2*x(1) + 18*sin(18*x(1));
        g(2,1) = 2*x(2) + 18*sin(18*x(2));
    case 4 % GENO The Economic Dispatch Problem
        g = zeros(4,1);
        g(1,1) = 3+3*0.000001*x(1)^2;
        g(2,1) = 2-0.000002/3*3/(x(2)^4);
    case 5 % GENO A Pressure Vessel Design Problem
        g = zeros(4,1);
        g(1,1) = 0.6224*0.0625*x(3)*x(4) + 3.1661*0.0625^2*x(4)*2*x(1) + 19.84*0.0625^2*x(3)*2*x(1);
        g(2,1) = 0.0625*1.7781*(x(3)^2);
        g(3,1) = 0.6224*x(1)*0.0625*x(4) + 1.7781*0.0625*x(2)*2*x(3) + 19.84*0.0625^2*(x(1)^2);
        g(4,1) = 0.6224*x(1)*x(3)*0.0625 + 3.1661*(x(1)^2)*0.0625^2;
    case 6 % GENO The Alkylation Process
        g = zeros(10,1);
        g(1,1) = 5.04;
        g(2,1) = 0.035;
        g(3,1) = 10;
        g(4,1) = -0.063*x(7);
        g(5,1) = 3.36;
        g(7,1) = -0.063*x(4);
    case 7 % GENO Decentralised Economic Planning
        g = zeros(7,1);
        g(1,1) = 2*(x(1)-10);
        g(2,1) = 10*(x(2)-12);
        g(3,1) = 4*x(3)^3;
        g(4,1) = 6*(x(4)-11);
        g(5,1) = 60*x(5)^5;
        g(6,1) = 14*x(6) -4*x(7)-10;
        g(7,1) = 4*x(7)^3-4*x(6)-8;
    case 8 % GENO Heat Exchanger Optimisation
        g = zeros(8,1);
        g(1,1) = 1;
        g(2,1) = 1;
        g(3,1) = 1;
    case 9 % GENO The Harvest Problem
        g = zeros(9,1);
        g(1,1) = -0.5*x(1)^(-0.5);
        g(2,1) = -0.5*x(2)^(-0.5);
        g(3,1) = -0.5*x(3)^(-0.5);
        g(4,1) = -0.5*x(4)^(-0.5);
    case 10 % GENO A Non-linear Resource Allocation Problem
        g = zeros(11,1);
        for j=1:5
            g(j,1) = 1;
            for i=1:5
                if j~= i
                    g(j,1) = g(j,1)*(1+i*x(i));
                else
                    g(j,1) = g(j,1)*i;
                end
            end
        end
        g = -g;
    case 11 % GENO Oligopolist Market Equilibrium Problem
        ci = [10;8;6;4;2];
        alphai = [1/1.2; 1/1.1; 1.00; 1/0.9; 1/0.8];
        beta = 1/1.1;
        Q = sum(x);
        pQ = (5000/Q).^(beta);
        g = -diag( ci + x.^alphai );
        g = g + pQ;
        g(1,:) = g(1,:) + x(1)*pQ;
        g(2,:) = g(2,:) + x(1)*pQ;
        g(3,:) = g(3,:) + x(1)*pQ;
        g(4,:) = g(4,:) + x(1)*pQ;
        g(5,:) = g(5,:) + x(1)*pQ;
    case 12 % GENO The Euclidean Compromise Solution I
        g = [ 2*x(1) ; 2*(x(1)-2) ];
    case 13 % GENO The Euclidean Compromise Solution II
        g = [...
            2*(x(1)-1),2*(x(2)-3); ...
            2*(x(1)-4),2*(x(2)-2) ];
    case 14 % GENO Efficient Portfolio Selection
        g = zeros(11,1); % Fix later
    case 15 % GENO A Dynamic Non-cooperative Game 1
        % u1(1) ... u1(T) u2(1) ... u2(T) .. x1(1) x1(T+1) x2(1) x2(T+1)
        g = zeros(22,1);
        T = 5;
        N = 2;
        g(T*N+T+1,1) = -1;
        g(T+1:2*T,1) = 1/T*x(T+1:2*T);
    case 16 % GENO A Dynamic Non-cooperative Game 2
        T = 5;
        N = 2;
        g = zeros(22,1);
        g(T*N+T+1,1) = g(T*N+T+1,1)+sign(x(T*N+T+1));
        g(end,1) = g(end,1)+sign(x(end));
        g(N*T+1:N*T+T,1) = g(N*T+1:N*T+T,1)+sign(x(N*T+1:N*T+T));
        g(N*T+2+T:end-1,1) = g(N*T+2+T:end-1,1)+sign(x(N*T+2+T:end-1));
        g(1:T,1) = g(1:T,1)+0.5*sign(x(1:T));
        g(T+1:T*N,1) = g(T+1:T*N,1)+0.5*sign(x(T+1:T*N));
    case 17 % GENO A Dynamic Non-cooperative Game 3
        s = 0.5;
        q = 0.95;
        T = 5;
        v = 0.5;
        g = zeros(11,1);
        g(end,1) = s*((1-v)*x(end)^(-v))*q^T;
        g(1:T,1) = q.^(0:T-1)'.*(1-v).*x(1:T).^(-v);
        g = -g;
    case 18 % GENO A Dynamic Non-cooperative Game 4
        q = 1;
        s = 1;
        r = 1;
        T = 5;
        g = zeros(11,1);
        g(end,1) = 2*q*x(end);
        g(T+1:end-1,1) = 2*s*x(T+1:end-1);
        g(1:T,1) = 2*r*x(1:T);
    case 19 % GENO A Dynamic Non-cooperative Game 5
        T = 5;
        g = zeros(11,1);
        for j=1:T
            g(j,1) = 1;
            for i=1:T
                if i~=j
                    g(j,1) = g(j,1)*(1+i*x(i));
                else
                    g(j,1) = g(j,1)*i;
                end
            end
        end
        g = -g;
    case 20 % GENO A Dynamic Non-cooperative Game 6
        T = 5;
        g = zeros(11,1);
        g(end) = sign(x(end)-1);
        g(1:T,1) = 3.^(-x(T+1:end-1));
        g(T+1:end-1,1) = -(2+x(1:T)).*3.^(-x(T+1:end-1))*log(3);
    case 21 % GENO A Dynamic Non-cooperative Game 7
        T = 5;
        g = zeros(11,1);
        g(1:T,1) = -0.5*x(1:T).^(-0.5);
end

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 080603  med  Switched to *Assign, cleaned