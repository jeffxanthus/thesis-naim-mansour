function [ x ] = BasisPursuit( A, b, Bcon, Ccon, thetaClipPos, thetaClipNeg )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Solves ?||x||1 + 1 (||Ax ? b||2)2 ?2
% where ? = ?	2log(p) % Basis pursuit with denoising
[m,n] = size(A); 
gamma = 0.003;
% Our QP is %Min??n (ui+vi)+1pTp
% subject to A(u ? v) ? Ip = b, u ? 0andv ? 0
% We split p = (p1 ? p2) in the actual formulation

H = sparse(2*m+2*n,2*m+2*n); 
H(1:2*m,1:2*m) = [eye(m) -eye(m); -eye(m) eye(m)];

%H2 = sparse(m+n,m+n);
%H2(1:m-1,1:m-1) = [eye(floor(m/2)) -eye(floor(m/2)); -eye(floor(m/2)) eye(floor(m/2))];

f = [zeros(2*m,1); gamma*ones(2*n,1)];
%f2 = [zeros(m,1); gamma*ones(n,1)];
% Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A; zeros(2*m,n) zeros(2*m,n) -eye(2*m,2*n)]; 
% beq = [b; -b; zeros(n,1); zeros(n,1)];

% Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A; zeros(2*m,n) zeros(2*m,n) -eye(2*m,2*n); zeros(m) zeros(m) Bcon -Bcon; zeros(m) zeros(m) -Ccon Ccon]; 
% beq = [b; -b; zeros(n,1); zeros(n,1); thetaClipPos; thetaClipNeg];

Aeq = [-speye(m) speye(m) A -A; speye(m) -speye(m) -A A; zeros(2*m,n) zeros(2*m,n) -speye(2*m,2*n); zeros(m) zeros(m) Bcon -Bcon; zeros(m) zeros(m) -Ccon Ccon]; 
%Aeq2 = [-speye(m,n) A; speye(m,n) -A; zeros(m,n) speye(m,n); zeros(m) Bcon; zeros(m) -Ccon];

beq = [b; -b; zeros(n,1); zeros(n,1); thetaClipPos; thetaClipNeg];
%beq2 = [b; -b; zeros(n,1); thetaClipPos; thetaClipNeg;];
%size(beq2)


% Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A]; 
% beq = [b; -b];


%A = [zeros(m) zeros(m) -eye(m,2*n)];
%B = zeros(m,1);

%sol = sqsolqp(H,Aeq,beq,f,1e-2);

%opt = optimset('TolPCG', 1e-2, 'TolX',1e-2, 'TolFun',1e-2, 'Algorithm', 'interior-point-convex');
%sol = quadprog(H,f,A,B,Aeq,beq,-inf,inf, [], opt); 
% x is (u-v)

x_L = ones(size(Aeq,2),1)*(-2000);
x_U = ones(size(Aeq,2),1)*1000;
b_U = ones(size(Aeq,1),1)*2000;

tic
[sol] = cplex(f, Aeq, x_L, x_U, beq, b_U, [], [], [], [], [],[], [], [], [], [], H);
toc

     
x = sol(2*m+1:2*m+n) - sol(2*m+n+1:2*m+2*n);

end

